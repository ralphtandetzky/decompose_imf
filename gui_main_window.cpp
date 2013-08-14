#include "gui_main_window.h"
#include "ui_gui_main_window.h"

#include "calculations.h"
#include "std_make_unique.h"
#include "optimize.h"
#include "sqr.h"

#include <complex>
#include <vector>
#include <array>
#include <iostream>

namespace gui {

MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent)
    , ui( std::make_unique<Ui::MainWindow>() )
{
    ui->setupUi(this);

    std::vector<double> f;
    const int nSamples = 15;

    for ( auto i = 0; i < nSamples; ++i )
    {
        f.push_back( sin(i*2*3.141592/nSamples*2) );
    }

    std::vector<std::vector<double>> swarm(300,
        std::vector<double>((f.size()+1)*2));
    std::minstd_rand rng;
    std::normal_distribution<> normal_dist;
    for ( auto & x : swarm )
        for ( auto & t : x )
            t = normal_dist(rng);

    auto display = []( const std::vector<double> & v )
    {
        for ( const auto & elem : v )
            printf( "%5d;", int(std::round(100*elem)));
        std::cout << std::endl;

    };

    int nIters = 10000;
    swarm = differentialEvolution(
        std::move(swarm),
        /*CO =*/ 0.1,
        /*DW =*/ 0.6,
        [&f]( const std::vector<double> & v ) -> double
        {
            return costFunction( f, v );
        },
        [&]( const decltype(swarm) & )-> bool
        {
            return --nIters <= 0;
        },
        [&]( const std::vector<double> & v, double cost )
        {
            std::cout << nIters << ' ' << cost << ' ';
            display(v);
        },
        rng );

    display(f);
    for ( const auto & row : swarm )
        display( calculateImfFromPairsOfReals(row) );
    std::cout << std::endl;
    for ( const auto & row : swarm )
        display( row );
}

MainWindow::~MainWindow()
{
}

} // namespace gui
