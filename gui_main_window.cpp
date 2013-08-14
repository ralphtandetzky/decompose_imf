#include "gui_main_window.h"
#include "ui_gui_main_window.h"

#include "calculations.h"
#include "cpp_utils/std_make_unique.h"
#include "optimize.h"
#include "cpp_utils/sqr.h"

#include <complex>
#include <vector>
#include <array>
#include <iostream>

namespace gui {

struct MainWindow::Impl
{
    Ui::MainWindow ui;
};

MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent)
    , m( std::make_unique<Impl>() )
{
    m->ui.setupUi(this);

}

MainWindow::~MainWindow()
{
}

void MainWindow::optimize()
{
    std::vector<double> f;
    const int nSamples = m->ui.nSamplesSpinBox->value();
    const int swarmSize = m->ui.swarmSizeSpinBox->value();
    const double crossOverProb = m->ui.coSpinBox->value();
    const double diffWeight    = m->ui.dwSpinBox->value();

    for ( auto i = 0; i < nSamples; ++i )
    {
        f.push_back( cos(i*2*3.141592/nSamples*2) );
    }

    std::vector<std::vector<double>> swarm( swarmSize,
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
        crossOverProb,
        diffWeight,
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


} // namespace gui
