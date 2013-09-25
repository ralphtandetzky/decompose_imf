#include "gui_main_window.h"
#include "ui_gui_main_window.h"

#include "calculations.h"
#include "cpp_utils/std_make_unique.h"
#include "cpp_utils/optimize.h"
#include "cpp_utils/sqr.h"
#include "cpp_utils/parallel_executor.h"
#include "cpp_utils/math_constants.h"

#include "qt_utils/serialize_props.h"

#include <complex>
#include <vector>
#include <array>
#include <iostream>
#include <thread>
#include <fstream>

#include <QtGui>


#include "cpp_utils/user_parameter.h"

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
    std::atomic<bool> cancelled;
    // single-threaded for background ops.
    cu::ParallelExecutor worker;
};


MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent)
    , m( std::make_unique<Impl>() )
{
    cu::RealUserParameter rp(0.1,std::string(),std::string(),std::string(),0.0,1.0,0.1,2,std::string(),1);

    m->ui.setupUi(this);
//    qu::createPropertySerializers( this->findChildren<QCheckBox*>(),
//                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QDoubleSpinBox*>(),
                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QComboBox*>(),
                                   std::back_inserter( m->serializers ) );
    qu::createPropertySerializers( this->findChildren<QSpinBox*>(),
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
    cancel();
    const int    nSamples      = m->ui.nSamplesSpinBox ->value();
    const int    swarmSize     = m->ui.swarmSizeSpinBox->value();
    const double crossOverProb = m->ui.coSpinBox       ->value();
    const double diffWeight    = m->ui.dwSpinBox       ->value();
    const double angleDev      = m->ui.aSpinBox        ->value()/360;
    const double amplitudeDev  = m->ui.bSpinBox        ->value();

    m->worker.addTask( [=]()
    {
        using cu::pi;
        m->cancelled = false;

        std::vector<double> f;

        for ( auto i = 0; i < nSamples; ++i )
        {
            f.push_back( cos(i*2*pi/nSamples*2) );
        }

        std::vector<std::vector<double>> swarm( swarmSize,
            std::vector<double>((f.size()+1)*2));
        std::minstd_rand rng;
        std::normal_distribution<> normal_dist;
        std::uniform_real_distribution<double> uniform(-pi,pi);
        for ( auto & x : swarm )
        {
            for ( size_t i = 0; i < x.size(); i+=2 )
            {
                x[i  ] = amplitudeDev*normal_dist(rng);
                x[i+1] = angleDev*uniform(rng) + i*pi/nSamples*2 + pi/2;
            }
        }

        auto display = []( const std::vector<double> & v )
        {
            for ( const auto & elem : v )
                printf( "%5d;", int(std::round(100*elem)));
            std::cout << std::endl;
        };

        int nIter = 0;

        const auto cost = [&f]( const std::vector<double> & v ) -> double
        {
            for ( size_t i = 1; i < v.size(); i+=2 )
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
            const int scale = 20;
            const int psize = v.size()*scale/2;
            QPixmap pixmap(psize,psize);
            {
                QPainter painter{&pixmap};
                painter.fillRect(0,0,psize,psize,Qt::black);
                QPolygonF reals, imags;
                for ( size_t i = 0; i < v.size(); i+=2 )
                {
                    reals << QPointF( scale*i/2, -scale*v[i  ] );
                    imags << QPointF( scale*i/2, -scale*std::remainder(v[i+1],2*pi) );
                }
                const auto imf = calculateImfFromPairsOfReals( v );
                QPolygonF fPoly, imfPoly;
                for ( size_t i = 0; i < f.size(); ++i )
                {
                    fPoly   << QPointF( scale*i, -scale*f  [i] );
                    imfPoly << QPointF( scale*i, -scale*imf[i] );
                }
                reals  .translate(0,psize/2);
                imags  .translate(0,psize/2);
                fPoly  .translate(0,psize/2);
                imfPoly.translate(0,psize/2);
                painter.setRenderHint( QPainter::Antialiasing );
                painter.setPen( Qt::darkGray );
                painter.drawLine( 0, psize/2, psize, psize/2);
                painter.drawLine( 0, psize/2-pi*scale, psize, psize/2-pi*scale);
                painter.drawLine( 0, psize/2+pi*scale, psize, psize/2+pi*scale);
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

        display(f);
        for ( const auto & row : swarm )
            display( calculateImfFromPairsOfReals(row) );
        std::cout << std::endl;
        for ( const auto & row : swarm )
            display( row );
    } );
}


void MainWindow::cancel()
{
    m->cancelled = true;
}


} // namespace gui
