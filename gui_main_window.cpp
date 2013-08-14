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
    const int nSamples = 10;

    for ( auto i = 0; i < nSamples; ++i )
    {
        f.push_back( sin(i*2*3.141592/nSamples) );
    }

    std::vector<std::vector<double>> swarm(30,std::vector<double>((f.size()+1)*2));
    std::minstd_rand rng;
    std::normal_distribution<> normal_dist;
    for ( auto & x : swarm )
        for ( auto & t : x )
            t = normal_dist(rng);

    int nIters = 1000;
    swarm = differentialEvolution(
        std::move(swarm),
        /*CO =*/ 0.5,
        /*DW =*/ 0.6,
        [&f]( const std::vector<double> & v ) -> double
        {
            return costFunction( f, v );
        },
        [&]( const decltype(swarm) & )-> bool
        {
            return --nIters <= 0;
        },
        rng );

    auto display = []( const std::vector<double> & v )
    {
        for ( const auto & elem : v )
            printf( "%5d;", int(std::round(100*elem)));
        std::cout << std::endl;

    };
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
