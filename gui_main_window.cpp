#include "gui_main_window.h"
#include "ui_gui_main_window.h"

#include "calculations.h"

#include "cpp_utils/formula_parser.h"
#include "cpp_utils/math_constants.h"
#include "cpp_utils/optimize.h"
#include "cpp_utils/parallel_executor.h"
#include "cpp_utils/sqr.h"
#include "cpp_utils/std_make_unique.h"
#include "cpp_utils/user_parameter.h"

#include "qt_utils/event_filter.h"
#include "qt_utils/serialize_props.h"

#include <array>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include <QtGui>
#include <QMessageBox>
#include <QEvent>
#include <QFileDialog>

#include <opencv/cv.h>

namespace gui {

struct MainWindow::Impl
{
    Impl()
        : worker{ 1 } // create onle one thread.
    {
    }

    ~Impl()
    {
    }

    /////////////////////
    // Gui thread data //
    /////////////////////

    // Contains Qt user interface elements.
    Ui::MainWindow ui;
    // Objects which help to load the values in the gui input widgets
    // during construction and to store them during destruction.
    std::vector<std::unique_ptr<qu::PropertySerializer>> serializers;
    // A table of different methods which find an approximate solution
    // to the optimization problem. The user can select the method in
    // the gui.
    std::map<std::string,std::function<
        std::vector<std::complex<double>>(
            const std::vector<double> &)>> initializers;
    // A vector of samples read from a file.
    std::vector<double> samples;


    ///////////////////////////////////////////////
    // Data shared between gui thread and worker //
    ///////////////////////////////////////////////

	struct SharedData
	{	
	    // Holds whether the currently running optimization task (if any)
	   	// has been cancelled.
	   	bool cancelled = true;
	    // Holds whether the currently running optimization task (if any)
	    // should proceed with calculating the imf of the current residue
	    // function.
	    bool shall_calculate_next_imf = false;
	};
	cu::Monitor<SharedData> shared;

    ///////////////////
    // Worker thread //
    ///////////////////

    // single-threaded for background ops.
    cu::ParallelExecutor worker;
};


MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent)
    , m( std::make_unique<Impl>() )
{
    m->ui.setupUi(this);

    // create table of functions to find initial approximations
    m->initializers["Zero"] =
        []( const std::vector<double> & f )
        { return std::vector<std::complex<double>>(f.size()+1); };
    m->initializers["Interpolate zeros"] =
        &getInitialApproximationByInterpolatingZeros;

    // display these methods in the gui.
    for ( const auto & initializer : m->initializers )
        m->ui.initApproxMethComboBox->addItem(
            QString::fromStdString(initializer.first) );

    // set up serializers
//    qu::createPropertySerializers( this->findChildren<QCheckBox*>(),
//                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QDoubleSpinBox*>(),
                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QComboBox*>(),
                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QSpinBox*>(),
                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QLineEdit*>(),
                                   std::back_inserter( m->serializers ) );

    // load serialized input widget entries from a settings file.
    std::ifstream file( "settings.txt" );
    readProperties( file, m->serializers );
}


MainWindow::~MainWindow()
{
    // stop any running optimization thread
    cancel();

    // store current values from gui input widget entries.
    std::ofstream file( "settings.txt" );
    writeProperties( file, m->serializers );
}


void MainWindow::optimize()
{
    // cancel the currently running operation (if any)
    cancel();

    // read values from gui
    const auto functionString = m->ui.functionLineEdit->text().toStdString();
    const auto xmin           = m->ui.xminSpinBox     ->value();
    const auto xmax           = m->ui.xmaxSpinBox     ->value();
    const auto nSamplesExpr   = m->ui.nSamplesSpinBox ->value();
    const auto swarmSize      = m->ui.swarmSizeSpinBox->value();
    const auto angleDevDegs   = m->ui.angleDevSpinBox ->value();
    const auto amplitudeDev   = m->ui.amplDevSpinBox  ->value();
    const auto crossOverProb  = m->ui.coSpinBox       ->value();
    const auto diffWeight     = m->ui.dwSpinBox       ->value();
    const auto initializer    = m->initializers.at(
                m->ui.initApproxMethComboBox->currentText().toStdString());
    const size_t nParams      = m->ui.nBaseFuncsSpinBox->value();
    const auto initSigmaUnits = m->ui.initSigmaSpinBox ->value();
    const auto initTauUnits   = m->ui.initTauSpinBox   ->value();
    const auto nodeDevUnits   = m->ui.nodeDevSpinBox   ->value();
    const auto sigmaDevUnits  = m->ui.sigmaDevSpinBox  ->value();
    const auto tauDevUnits    = m->ui.tauDevSpinBox    ->value();


    // build expression tree for target function
    const auto expression = std::make_shared<cu::ExpressionTree>();
    const auto nParsedChars = expression->parse( functionString );
    if ( nParsedChars < functionString.size() )
    {
        std::stringstream ss;
        ss << "Target function is not valid. "
              "There seems to be an error here: \n";
        ss << std::string( begin(functionString),
                           begin(functionString)+nParsedChars );
        ss << std::endl << ">>\t";
        ss << std::string( begin(functionString)+nParsedChars,
                           end(functionString) );
        QMessageBox msgBox;

        msgBox.setText( QString::fromStdString(ss.str()) );
        msgBox.exec();
        return;
    }

    // calculate the target function from the expression
    auto f = std::vector<double>{};
    const auto tab = m->ui.samplesTabWidget->currentWidget();
    if ( tab == m->ui.fromExpressionTab )
    {
        for ( auto i = 0; i < nSamplesExpr; ++i )
            f.push_back( xmin + (xmax-xmin)*i/(nSamplesExpr-1) );
        f = expression->evaluate( f );
    }
    else if ( tab == m->ui.fromFileTab )
    {
        f = m->samples;
    }
    else
        CU_THROW( "No tab is open for selecting the samples." );
    const auto nSamples = f.size();

    // concurrently launch optimization on the worker thread
    m->worker.addTask( [=]() mutable // only f is mutated
    {
        using cu::pi;
        m->shared( []( Impl::SharedData & shared )
        { 
        	shared.cancelled = false; 
        	shared.shall_calculate_next_imf = false; 
        } );

        auto done = false;

        while ( !done )
        {
            // calculate an initial approximation and swarm
            const auto initApprox = initializer(f);

            // calculate base with equidistant center points of logistic functions
            auto logisticBase = std::vector<std::vector<double>>{};
            auto nodes = std::vector<double>{};
            const auto factor = nSamples / (xmax-xmin);
            const auto initSigma = initSigmaUnits * factor;
            for ( auto i = size_t{0}; i < nParams; ++i )
            {
                nodes.push_back( (i+.5)*initApprox.size()/nParams );
                logisticBase.push_back( getSamplesFromLogisticFunctionBase(
                    { 1., nodes.back() }, initSigma, initApprox.size() ) );
            }
            auto logisticBaseMat = cv::Mat(
                        logisticBase.front().size(), logisticBase.size(),
                        CV_64FC1, cv::Scalar::all(0) );
            for ( auto row = 0; row < logisticBaseMat.rows; ++row )
            {
                for ( auto col = 0; col < logisticBaseMat.cols; ++col )
                {
                    logisticBaseMat.at<double>(row,col) =
                            logisticBase.at(col).at(row);
                }
            }

            // calculate the element closest to initApprox that lies in the
            // space spanned by the base elements. Also calculate the
            // coefficients of that minimum element with respect to the given
            // base
            auto initApproxMat = cv::Mat(
                        initApprox.size(), 1, CV_64FC1, cv::Scalar::all(0) );
            for ( auto row = 0; row < initApproxMat.rows; ++row )
            {
                initApproxMat.at<double>( row ) = initApprox.at(row).imag();
            }
            const auto invLogisticBase =
                    cv::Mat{ logisticBaseMat.inv( cv::DECOMP_SVD ) };
            const auto bestApproxMat =
                    cv::Mat{ invLogisticBase * initApproxMat };

            auto swarm = std::vector<std::vector<double>>( swarmSize );
            auto rng = std::mt19937{};
            {
                auto normal_dist = std::normal_distribution<>{};
                auto uniform = std::uniform_real_distribution<double>{-1,1};
                for ( auto & x : swarm )
                {
                    const auto sigmaDev  = sigmaDevUnits  * factor;
                    const auto tauDev    = tauDevUnits    * factor;
                    const auto initTau   = initTauUnits   * factor;
                    const auto nodeDev   = nodeDevUnits   * factor;
                    const auto angleDev  = angleDevDegs/180*cu::pi;
                    for ( auto i = size_t{0}; i < nodes.size(); ++i )
                    {
                        x.push_back( amplitudeDev*normal_dist(rng) );
                        x.push_back( nodes[i] + nodeDev*normal_dist(rng) );
                    }
                    for ( auto i = size_t{0}; i < nodes.size(); ++i )
                    {
                        x.push_back( bestApproxMat.at<double>( i ) +
                                     angleDev*uniform(rng) );
                        x.push_back( nodes[i] );
                    }
                    x.push_back( initSigma + sigmaDev*normal_dist(rng) );
                    x.push_back( initTau + tauDev*normal_dist(rng) );
                }
            }

            // cost function for optimization
            const auto cost = [&f, nSamples]( std::vector<double> v ) -> double
            {
                return costFunction( f,
                    getSamplesFromParams( std::move(v), nSamples ) );
            };

            // This variable is shared between 'shallTerminate' and 'sendBestFit'.
            auto nIter = 0;

            // function which returns whether the
            // optimization algorithm shall terminate.
            const auto shallTerminate = [&]( const decltype(swarm) & )-> bool
            {
                ++nIter;
                return m->shared( [&]( Impl::SharedData & shared ) -> bool
                {
                	if ( shared.shall_calculate_next_imf )
                	{
                		shared.shall_calculate_next_imf = false;
                		nIter = 0;
                		return true;
                	}
                	if ( shared.cancelled )
                	{
                		done = true;
                		return true;
                	}
                	return false;
                } );
            };

            // function which is called by the optimization algorithm
            // each time the best fit is improved.
            auto bestParams = std::vector<double>{};
            const auto sendBestFit = [&](
                    const std::vector<double> & v_, double cost )
            {
                bestParams = v_;
                const auto v = getSamplesFromParams( v_, nSamples );

                // console output
                std::cout << nIter << ' ' << cost << ' ' << std::endl;
                //for ( const auto & elem : v )
                //    printf( "%5d;", int(std::round(100*elem)));
                //std::cout << std::endl;

                const auto psize = 600.;
                const auto xscale = psize*2/(v.size()-2);
                const auto yscale = 20;
                QImage image{(int)psize,(int)psize,QImage::Format_RGB32};
                {
                    QPainter painter{&image};
                    painter.fillRect(0,0,psize,psize,Qt::black);
                    auto reals = QPolygonF{};
                    auto imags = QPolygonF{};
                    for ( auto i = size_t{0}; i < v.size(); i+=2 )
                    {
                        reals << QPointF( xscale*i/2, -yscale*v[i  ] );
                        imags << QPointF( xscale*i/2, -yscale*std::remainder(v[i+1],2*pi) );
                    }
                    const auto imf = calculateImfFromPairsOfReals( v );
                    auto fPoly   = QPolygonF{};
                    auto imfPoly = QPolygonF{};
                    for ( auto i = size_t{0}; i < f.size(); ++i )
                    {
                        fPoly   << QPointF( xscale*i, -yscale*f  [i] );
                        imfPoly << QPointF( xscale*i, -yscale*imf[i] );
                    }
                    reals  .translate(0,psize/2);
                    imags  .translate(0,psize/2);
                    fPoly  .translate(0,psize/2);
                    imfPoly.translate(0,psize/2);
                    painter.setRenderHint( QPainter::Antialiasing );
                    painter.setPen( Qt::darkGray );
                    painter.drawLine( 0, psize/2, psize, psize/2);
                    painter.drawLine( 0, psize/2-pi*yscale, psize, psize/2-pi*yscale);
                    painter.drawLine( 0, psize/2+pi*yscale, psize, psize/2+pi*yscale);
                    painter.setPen( Qt::green );
                    painter.drawPolyline( reals );
                    painter.setPen( Qt::magenta );
                    painter.drawPolyline( imags );
                    painter.setPen( Qt::white );
                    painter.drawPolyline( fPoly );
                    painter.setPen( Qt::yellow );
                    painter.drawPolyline( imfPoly );
                }
                QMetaObject::invokeMethod( this, "setDisplay",
                                           Q_ARG(QImage,image));
            };

            // perform the optimization.
            swarm = cu::differentialEvolution(
                std::move(swarm), crossOverProb, diffWeight,
                cost, shallTerminate, sendBestFit, rng );

            const auto bestImf = calculateImfFromPairsOfReals(
                        getSamplesFromParams( bestParams, nSamples ) );
            cu::for_each( begin(f), end(f), bestImf.begin(), bestImf.end(),
                          []( double & lhs, double rhs ){ lhs -= rhs; } );
        } // while loop
    } ); // task for worker
}


void MainWindow::cancel()
{
	m->shared( []( Impl::SharedData & shared )
	{ 
		shared.cancelled = true; 
	} );
}


void MainWindow::calculateNextImf()
{
	m->shared( []( Impl::SharedData & shared )
	{ 
		shared.shall_calculate_next_imf = true; 
	} );
}


void MainWindow::selectSamplesFile()
try
{
    const auto lineEdit = m->ui.samplesFileLineEdit;
    const auto qFileName = QFileDialog::getOpenFileName();
    const auto fileName = qFileName.toStdString();
    std::ifstream file{ fileName };
    if ( !file )
        CU_THROW( "Could not open the file \"" + fileName + "\"." );
    auto vals = std::vector<double>( std::istream_iterator<double>(file),
                                     std::istream_iterator<double>() );
    if ( file.bad() )
        CU_THROW( "The file \"" + fileName +
                  "\" could not be read." );
    if ( vals.empty() )
        CU_THROW( "The file \"" + fileName +
                  "\" does not contain samples." );
    if ( !file.eof() )
        CU_THROW( "The end of the file \"" + fileName +
                  "\" has not been reached." );

    lineEdit->setText( qFileName );
    m->samples = std::move(vals);
}
catch (...)
{
    CU_THROW( "Could not read samples file successfully." );
}


void MainWindow::setDisplay( const QImage & image )
{
    m->ui.graphDisplay->setPixmap( QPixmap::fromImage(image) );
}

} // namespace gui
