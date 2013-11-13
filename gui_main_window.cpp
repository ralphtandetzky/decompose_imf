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

namespace gui {

struct MainWindow::Impl
{
    Impl()
        : cancelled(true)
        , worker( 1 )
    {
    }

    ~Impl()
    {
    }

    Ui::MainWindow ui;
    std::vector<std::unique_ptr<qu::PropertySerializer>> serializers;
    std::map<std::string,std::function<
        std::vector<std::complex<double>>(
            const std::vector<double> &)>> initializers;
    std::atomic<bool> cancelled;
    // single-threaded for background ops.
    cu::ParallelExecutor worker;
};


MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent)
    , m( std::make_unique<Impl>() )
{
    m->initializers["Zero"] =
        []( const std::vector<double> & f )
        { return std::vector<std::complex<double>>(f.size()+1); };
    m->initializers["Interpolate zeros"] =
        &getInitialApproximationByInterpolatingZeros;
    m->ui.setupUi(this);

    for ( const auto & initializer : m->initializers )
        m->ui.initApproxMethComboBox->addItem(
            QString::fromStdString(initializer.first) );

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

    std::ifstream file( "settings.txt" );
    readProperties( file, m->serializers );
}


MainWindow::~MainWindow()
{
    cancel();
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
    const auto nSamples       = m->ui.nSamplesSpinBox ->value();
    const auto swarmSize      = m->ui.swarmSizeSpinBox->value();
    const auto angleDev       = m->ui.angleDevSpinBox ->value()/180*cu::pi;
    const auto amplitudeDev   = m->ui.amplDevSpinBox  ->value();
    const auto crossOverProb  = m->ui.coSpinBox       ->value();
    const auto diffWeight     = m->ui.dwSpinBox       ->value();
    const auto initializer    = m->initializers.at(
                m->ui.initApproxMethComboBox->currentText().toStdString());

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

    // concurrently launch optimization on the worker thread
    m->worker.addTask( [=]()
    {
        using cu::pi;
        m->cancelled = false;

        auto f = std::vector<double>{};

        for ( auto i = 0; i < nSamples; ++i )
            f.push_back( xmin + (xmax-xmin)*i/(nSamples-1) );
        f = expression->evaluate( f );

        const auto initApprox = initializer(f);
        auto swarm = std::vector<std::vector<double>>( swarmSize,
            std::vector<double>(2*initApprox.size()) );
        auto rng = std::mt19937{};
        {
            auto normal_dist = std::normal_distribution<>{};
            auto uniform = std::uniform_real_distribution<double>{-1,1};
            for ( auto & x : swarm )
            {
                for ( size_t i = 0; i < initApprox.size(); ++i )
                {
                    x[2*i  ] = initApprox[i].real()+amplitudeDev*normal_dist(rng);
                    x[2*i+1] = initApprox[i].imag()+angleDev*uniform(rng);
                }
            }
        }

        auto display = []( const std::vector<double> & v )
        {
            for ( const auto & elem : v )
                printf( "%5d;", int(std::round(100*elem)));
            std::cout << std::endl;
        };

        auto nIter = 0;

        const auto cost = [&f]( const std::vector<double> & v ) -> double
        {
            for ( auto i = size_t{1}; i < v.size(); i+=2 )
                const_cast<double &>(v[i]) = remainder( v[i], 2*pi );
            return costFunction( f, v );
        };

        const auto shallTerminate = [&]( const decltype(swarm) & )-> bool
        {
            ++nIter;
            return m->cancelled.load();
        };

        const auto sendBestFit = [&]( const std::vector<double> & v, double cost )
        {
            std::cout << nIter << ' ' << cost << ' ' << std::endl;
            display(v);
            const auto psize = 600.;
            const auto xscale = psize*2/(v.size()-2);
            const auto yscale = 20;
            QPixmap pixmap{(int)psize,(int)psize};
            {
                QPainter painter{&pixmap};
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
            QMetaObject::invokeMethod( m->ui.graphDisplay, "setPixmap",
                                       Q_ARG(QPixmap,pixmap));
        };

        swarm = cu::differentialEvolution(
            std::move(swarm), crossOverProb, diffWeight,
            cost, shallTerminate, sendBestFit, rng );

//        swarm = nelderMead( swarm, cost, shallTerminate, sendBestFit );

//        display(f);
//        for ( const auto & row : swarm )
//            display( calculateImfFromPairsOfReals(row) );
//        std::cout << std::endl;
//        for ( const auto & row : swarm )
//            display( row );
    } );
}


void MainWindow::cancel()
{
    m->cancelled = true;
}


} // namespace gui
