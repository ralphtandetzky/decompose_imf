#pragma once

#include <complex>
#include <vector>

std::vector<std::complex<double>> groupPairsToComplex(
        const std::vector<double> & pairsOfReals );

std::vector<std::complex<double>> calculateSigmaFunctionFromSigmaSequence(
        std::vector<std::complex<double>> sigma );

std::vector<double> calculateImfFromSigmaFunction(
        const std::vector<std::complex<double>> & sigma );

std::vector<double> calculateImfFromPairsOfReals(
        const std::vector<double> & input );

double sumOfSquaresOfDifference(
        const std::vector<double> & lhs,
        const std::vector<double> & rhs );

double costFunction( const std::vector<double> & f,
                     const std::vector<double> & pairsOfReals );
