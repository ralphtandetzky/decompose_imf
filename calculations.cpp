#include <cassert>
#include <complex>
#include <vector>

#include "sqr.h"


std::vector<std::complex<double>>
    groupPairsToComplex(
        const std::vector<double> & input )
{
    const auto size = input.size();
    assert( size % 2 == 0 );
    std::vector<std::complex<double>> result;
    for ( size_t i = 0; i < size; i+=2 )
        result.push_back( std::complex<double>(input[i],input[i+1]));

    return result;
}


std::vector<std::complex<double>>
    calculateSigmaFunctionFromSigmaSequence(
        std::vector<std::complex<double>> sigma )
{
    const auto size = sigma.size();
    for ( size_t i = 1; i < size; ++i )
        sigma[i-1]=(sigma[i-1]+sigma[i])/2.;
    sigma.pop_back();
    return sigma;
}


std::vector<double> calculateImfFromSigmaFunction(
        const std::vector<std::complex<double>> & sigma )
{
    std::vector<double> result;
    for ( const auto & z : sigma )
        result.push_back( exp(z.real())*cos(z.imag()));
    return result;
}


std::vector<double> calculateImfFromPairsOfReals(
        const std::vector<double> & input )
{
    return calculateImfFromSigmaFunction(
        calculateSigmaFunctionFromSigmaSequence(
            groupPairsToComplex( input ) ) );
}


double sumOfSquaresOfDifference(
        const std::vector<double> & lhs,
        const std::vector<double> & rhs )
{
    assert( lhs.size() == rhs.size() );
    const auto size = lhs.size();
    double sum = 0;
    for ( size_t i = 0; i < size; ++i )
        sum += sqr(lhs[i]-rhs[i]);
    return sum;
}


double costFunction( const std::vector<double> & f,
                     const std::vector<double> & pairsOfReals )
{
    return sumOfSquaresOfDifference( f,
        calculateImfFromPairsOfReals(pairsOfReals) );
}
