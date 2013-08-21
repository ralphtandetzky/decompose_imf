#include "cpp_utils/sqr.h"
#include "calculations.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <vector>


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
    const auto pi = 3.14159265358979;
    const auto tau = derive( sigma );
    const auto size = sigma.size();
    for ( size_t i = 1; i < size; ++i )
    {
        sigma[i-1]=(sigma[i-1]+sigma[i])/2.;
        if ( abs(tau[i-1].imag()) > pi )
            sigma[i-1].imag( sigma[i-1].imag()+pi );
    }
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
            groupPairsToComplex(input)));
}


double sumOfSquaresOfDifference(
        const std::vector<double> & lhs,
        const std::vector<double> & rhs )
{
    assert( lhs.size() == rhs.size() );
    const auto size = lhs.size();
    double sum = 0;
    for ( size_t i = 0; i < size; ++i )
        sum += cu::sqr(lhs[i]-rhs[i]);
    return sum;
}


double costFunction( const std::vector<double> & f,
                     const std::vector<double> & pairsOfReals )
{
    auto sigma_seq = groupPairsToComplex( pairsOfReals );

    return sumOfSquaresOfDifference( f,
        calculateImfFromSigmaFunction(
            calculateSigmaFunctionFromSigmaSequence(sigma_seq) ) )
        + boundaryCondition( sigma_seq );
}


std::vector<std::complex<double>>
    derive( std::vector<std::complex<double>> f )
{
    const auto size = f.size();
    for ( size_t i = 1; i < size; ++i )
        f[i-1] = f[i] - f[i-1];
    f.pop_back();
    return f;
}


double boundaryCondition( std::vector<std::complex<double>> sigma_seq )
{
    const double pi = 3.14159265358979;
    auto tau = derive( std::move(sigma_seq) );
    for_each( begin(sigma_seq), end(sigma_seq), [pi](std::complex<double>&c)
    { c.imag( std::remainder( c.imag(), 2*pi ) ); } );
    double result = 0;
    for ( size_t i = 1; i < tau.size(); ++i )
    {
        const auto lhs = abs(tau[i]-tau[i-1]);
        const auto rhs = cu::sqr( std::max(0.,
            std::min(tau[i].imag(),tau[i-1].imag()) ) );
        if ( lhs > rhs )
            result += lhs-rhs+1;
    }
    for ( const auto & t : tau )
    {
        if ( t.imag() <= 0 )
        {
            result += 1 - t.imag();
        }
    }
    return result;
}
