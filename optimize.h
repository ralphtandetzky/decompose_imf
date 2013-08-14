#pragma once

#include <cassert>
#include <random>

template < typename Container
         , typename CostFunction
         , typename ShallTerminateFunctor
         , typename SendBestFitFunctor
         , typename RandomNumberGenerator
           >
Container differentialEvolution(
        Container swarm,
        const double crossOverProbability,
        const double differentialWeight,
        CostFunction costFunction,
        ShallTerminateFunctor shallTerminate,
        SendBestFitFunctor sendBestFit,
        RandomNumberGenerator & rng )
{
    assert( swarm.size() >= 4 );
    const auto N = swarm[0].size();
    for ( auto & x : swarm )
        assert( x.size() == N );
    typedef decltype(costFunction(swarm[0])) Cost;
    std::vector<Cost> costs;
    std::transform( begin(swarm), end(swarm),
        back_inserter(costs), costFunction );
    auto lowestCostIndex = std::max_element(
        begin(swarm), end(swarm) ) - begin(swarm);
    sendBestFit( swarm[lowestCostIndex],
                 costs[lowestCostIndex] );

    typedef typename Container::size_type size_type;
    std::uniform_int_distribution<> disX(0, swarm.size()-1);
    std::uniform_int_distribution<size_type> disR(0, N-1);
    std::uniform_real_distribution<> uniform;

    while ( !shallTerminate( swarm ) )
    {
        for ( size_type x = 0; x < swarm.size(); ++x )
        {
            size_type a = 0, b = 0, c = 0;
            do a = disX(rng); while ( a == x );
            do b = disX(rng); while ( b == x || b == a );
            do c = disX(rng); while ( c == x || c == a || c == b );
            const auto R = disR(rng);
            auto y = swarm[x];
            for ( size_type i = 0; i < y.size(); ++i )
            {
                auto r = uniform(rng);
                if ( r < crossOverProbability || i == R )
                    y[i] = swarm[a][i] + differentialWeight *
                            ( swarm[b][i] - swarm[c][i] );
            }
            const auto cost = costFunction(y);
            if ( cost < costs[x] )
            {
                swarm[x] = y;
                costs[x] = cost;
                if ( cost < costs[lowestCostIndex] )
                {
                    lowestCostIndex = x;
                    sendBestFit( swarm[lowestCostIndex],
                                 costs[lowestCostIndex] );
                }
            }
        }
    }

    return swarm;
}
