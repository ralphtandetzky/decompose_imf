#include "processing.h"

#include "cpp_utils/extract_by_line.h"
#include "cpp_utils/more_algorithms.h"
#include "cpp_utils/exception.h"

#include <algorithm>
#include <iterator>
#include <sstream>


std::vector<double> processSamples(
        std::vector<double> samples
        , const std::string & instructions
        , const std::map<std::string,
            std::function<std::vector<double>(
                const std::vector<double> &,
                std::vector<double>)> > & functions
        )
{
    const auto processing =
        cu::extractByLine( std::istringstream( instructions ) );

    for ( auto line : processing )
    {
        try
        {
            line = cu::trim( line );
            if ( line.empty() )
                continue;
            std::istringstream is(line);
            std::string funcName;
            is >> funcName;
            if ( functions.count(funcName) == 0 )
            {
                const auto message =
                        "The function '" + funcName +
                        "' is unknown.";
                auto minDist = funcName.size() / 2 + 1;
                auto minDistName = std::string{};
                for ( const auto & x : functions )
                {
                    const auto dist = cu::levenshteinDistance(
                                funcName, x.first );
                    if ( dist < minDist )
                    {
                        minDist = dist;
                        minDistName = x.first;
                    }
                }
                // no appropriate match?
                if ( minDistName.empty() )
                    CU_THROW( message );
                CU_THROW( message + " Did you mean '" +
                          minDistName + "'?" );
            }
            std::vector<double> args;
            std::copy( std::istream_iterator<double>(is),
                       std::istream_iterator<double>(),
                       std::back_inserter(args) );
            if ( !is.eof() )
                CU_THROW( "Could not parse the arguments up until "
                          "the end of the line." );
            samples = functions.at(funcName)(args,samples);
        }
        catch (...)
        {
            CU_THROW( "The line '" + line + "' could not be executed." );
        }
    }

    return samples;
}


std::map<std::string,
    std::function<std::vector<double>(
        const std::vector<double> &,
        std::vector<double>
    )>> createProcessingFunctions()
{
    auto functions = std::map<std::string,
        std::function<std::vector<double>(
            const std::vector<double> &,
            std::vector<double>
        )>>{};

    functions["box_filter"] =
        []( const std::vector<double> & args, std::vector<double> samples )
    {
        if ( args.size() != 1 )
            CU_THROW( "The 'box_filter' preprocessing step expects "
                      "exactly one argument, not " +
                      std::to_string(args.size()) + "." );
        const auto width = size_t(args.front()+0.5);
        if ( width == 0 )
            CU_THROW( "Zero width is not a valid argument for the "
                      "box_filter preprocessing step." );
        const auto nSamples = samples.size();
        if ( width >= nSamples )
            CU_THROW( "Width " + std::to_string(width)+ " of box_filter "
                      "is too large for " + std::to_string(nSamples) +
                      " samples." );
        std::partial_sum( begin(samples), end(samples), begin(samples) );
        samples.insert( begin(samples), 0 );
        for ( auto i = size_t{0}; i + width < samples.size(); ++i )
            samples[i] = ( samples[i+width] - samples[i] ) / width;
        samples.resize( nSamples+1-width );
        return samples;
    };
    functions["clip"] =
        []( const std::vector<double> & args, std::vector<double> samples )
    {
        if ( args.size() != 2 )
            CU_THROW( "The 'clip' preprocessing step expects "
                      "exactly two arguments, not " +
                      std::to_string(args.size()) + "." );
        const auto first = size_t(args[0]+0.5);
        const auto last  = size_t(args[1]+0.5);
        if ( first >= last )
            CU_THROW( "The second argument of 'clip' must be "
                      "greater than the first one. "
                      "The first argument is " + std::to_string(first) +
                      " and the second is " + std::to_string(last) +
                      "." );
        if ( last > samples.size() )
            CU_THROW( "The upper bound " + std::to_string(last) +
                      " is greater than the number of samples " +
                      std::to_string(samples.size()) + "." );
        return std::vector<double>( samples.begin()+first,
                                    samples.begin()+last );
    };
    functions["high_pass"] =
        []( const std::vector<double> & args, std::vector<double> samples )
    {
        if ( args.size() != 1 )
            CU_THROW( "The 'high_pass' preprocessing step expects "
                      "exactly one argument, not " +
                      std::to_string(args.size()) + "." );
        const auto factor = args.front();
        if ( factor <= 0 )
            CU_THROW( "Invalid argument " + std::to_string(factor) +
                      " for 'high_pass'. Argument must be positive." );
        const auto rec = 1./(1.+factor);
        double x = samples.front();
        const auto lambda = [&]( double & s ) { s -= x = (s+factor*x)*rec; };
        std::for_each( samples. begin(), samples. end(), lambda );
        x = samples.back();
        std::for_each( samples.rbegin(), samples.rend(), lambda );
        return samples;
    };
    functions["low_pass"] =
        []( const std::vector<double> & args, std::vector<double> samples )
    {
        if ( args.size() != 1 )
            CU_THROW( "The 'low_pass' preprocessing step expects "
                      "exactly one argument, not " +
                      std::to_string(args.size()) + "." );
        const auto factor = args.front();
        if ( factor <= 0 )
            CU_THROW( "Invalid argument " + std::to_string(factor) +
                      " for 'low_pass'. Argument must be positive." );
        const auto rec = 1./(1.+factor);
        double x = samples.front();
        const auto lambda = [&]( double & s ) { s = x = (s+factor*x)*rec; };
        std::for_each( samples. begin(), samples. end(), lambda );
        std::for_each( samples.rbegin(), samples.rend(), lambda );
        return samples;
    };
    functions["mul"] =
        []( const std::vector<double> & args, std::vector<double> samples )
    {
        if ( args.size() != 1 )
            CU_THROW( "The 'mul' preprocessing step expects "
                      "exactly one argument, not " +
                      std::to_string(args.size()) + "." );
        const auto factor = args.front();
        std::for_each( begin(samples), end(samples),
                       [&]( double & s ){ s *= factor; } );
        return samples;
    };
    functions["zero_average"] =
        []( const std::vector<double> & args, std::vector<double> samples )
    {
        if ( args.size() != 0 )
            CU_THROW( "The 'zero_average' preprocessing step expects "
                      "no arguments, but " +
                      std::to_string(args.size()) + " have been passed." );
        const auto average =
                std::accumulate( begin(samples),
                                 end(samples), 0. ) / samples.size();
        std::for_each( begin(samples), end(samples),
                       [&]( double & s ){ s -= average; } );
        return samples;
    };
    return functions;
}