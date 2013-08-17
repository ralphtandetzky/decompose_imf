#include "gui_main_window.h"
#include "ui_gui_main_window.h"

#include "calculations.h"
#include "cpp_utils/std_make_unique.h"
#include "optimize.h"
#include "cpp_utils/sqr.h"
#include "cpp_utils/parallel_executor.h"

#include <complex>
#include <vector>
#include <array>
#include <iostream>
#include <thread>

#include <QtGui>

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
    std::atomic<bool> cancelled;
    // single-threaded for background ops.
    ParallelExecutor worker;
};

MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent)
    , m( std::make_unique<Impl>() )
{
    m->ui.setupUi(this);
}

MainWindow::~MainWindow()
{
    cancel();
}

void MainWindow::optimize()
{
    cancel();
    const int    nSamples      = m->ui.nSamplesSpinBox ->value();
    const int    swarmSize     = m->ui.swarmSizeSpinBox->value();
    const double crossOverProb = m->ui.coSpinBox       ->value();
    const double diffWeight    = m->ui.dwSpinBox       ->value();

    m->worker.addTask( [=]()
    {
        m->cancelled = false;

        std::vector<double> f;

        for ( auto i = 0; i < nSamples; ++i )
        {
            f.push_back( cos(i*2*3.141592/nSamples*2) );
        }

        std::vector<std::vector<double>> swarm( swarmSize,
            std::vector<double>((f.size()+1)*2));
        std::minstd_rand rng;
        std::normal_distribution<> normal_dist;
        for ( auto & x : swarm )
        {
            for ( size_t i = 0; i < x.size(); i+=2 )
            {
                x[i  ] = 1*normal_dist(rng);
                x[i+1] = 1*normal_dist(rng) + i*3.141592/nSamples*2;
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
            QPixmap pixmap(300,300);
            {
                QPainter painter{&pixmap};
                painter.fillRect(0,0,300,300,Qt::black);
                QPolygonF reals, imags;
                for ( size_t i = 0; i < v.size(); i+=2 )
                {
                    reals << QPointF( 5*i, -10*v[i  ] );
                    imags << QPointF( 5*i, -10*v[i+1] );
                }
                const auto imf = calculateImfFromPairsOfReals( v );
                QPolygonF fPoly, imfPoly;
                for ( size_t i = 0; i < f.size(); ++i )
                {
                    fPoly   << QPointF( 10*i, -10*f  [i] );
                    imfPoly << QPointF( 10*i, -10*imf[i] );
                }
                reals  .translate(0,150);
                imags  .translate(0,150);
                fPoly  .translate(0,150);
                imfPoly.translate(0,150);
                painter.setRenderHint( QPainter::Antialiasing );
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

        swarm = differentialEvolution(
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
