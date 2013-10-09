#include "calculations.h"

#include "cpp_utils/sqr.h"
#include "cpp_utils/math_constants.h"

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
        result.push_back( exp(z.real())*sin(z.imag()));
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
    for_each( begin(tau), end(tau), [pi](std::complex<double>&c)
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


std::vector<std::complex<double>> getInitialApproximationByInterpolatingZeros(
    const std::vector<double> & f )
{
    using cu::pi;

    std::vector<std::complex<double>> result( f.size()+1, 0 );
    if ( f.empty() )
    {
        return result;
    }

    // calculate zeros and extreme values between them
    std::vector<double> zeros;
    for ( size_t i = 1; i != f.size(); ++i )
    {
        // crossing zero?
        if ( ( f[i-1] < 0 && 0 < f[i  ] ) ||
             ( f[i  ] < 0 && 0 < f[i-1] ) )
            zeros.push_back( i - f[i]/(f[i]-f[i-1])+0.5 );
    }

    // no zeros?
    if ( zeros.empty() )
    {
        // assign imaginary parts
        result.assign( f.size()+1,
            std::complex<double>(0,f[0]<0?-pi/2:pi/2) );

        return result;
    }

    // assign imaginary parts
    result.front().imag( f.front() < 0 ? -pi/2 : pi/2 );
    for ( size_t i = 0; i < zeros.front(); ++i )
    {
        result[i].imag( result.front().imag() + i*(pi/2)/zeros.front() );
    }
    for ( size_t z = 1; z < zeros.size(); ++z )
    {
        for ( size_t i = zeros[z-1]; i < zeros[z]; ++i )
        {
            result[i].imag( result.front().imag()-pi/2
                            +pi*(z+(i-zeros[z-1])/(zeros[z]-zeros[z-1])) );
        }
    }
    for ( size_t i = zeros.back(); i < result.size(); ++i )
    {
        result[i].imag( result.front().imag()-pi/2
                        +pi*(zeros.size()+(i-zeros.back())/
                             (result.size()-zeros.back())) );
    }

    return result;
}
