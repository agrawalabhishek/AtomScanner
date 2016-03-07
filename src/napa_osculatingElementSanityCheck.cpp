/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

// This program will make use of the ATOM solver to construct a transfer trajectory
// between two points in space. Only a single transfer segment is simulated here, that is, no multitargeting.
// There is just one departure and one arrival object considered in the simulation. The object TLEs are
// taken from a catalog file. The user can specify which object will be departure and which will be arrival
// along with the time of flight. 

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <iterator>

#include <libsgp4/Globals.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <Atom/atom.hpp>
#include <Atom/convertCartesianStateToTwoLineElements.hpp>

#include <Astro/orbitalElementConversions.hpp>

#include <SML/sml.hpp>
#include <SML/constants.hpp>
#include <SML/basicFunctions.hpp>
#include <SML/linearAlgebra.hpp>

#include <boost/array.hpp>

#include </usr/local/abhi/pykep/src/lambert_problem.cpp>
#include </usr/local/abhi/pykep/src/lambert_problem.h>
#include </usr/local/abhi/pykep/src/keplerian_toolbox.h>


typedef double Real;
typedef std::vector < Real > Vector6;
typedef std::vector < Real > Vector3;
typedef std::vector < Real > Vector2;
typedef std::vector < std::vector < Real > > Vector2D;
typedef boost::array < Real, 3 > array3; 

//! Remove newline characters from string.
void removeNewline( std::string& string )
{
    string.erase( std::remove( string.begin( ), string.end( ), '\r' ), string.end( ) );
    string.erase( std::remove( string.begin( ), string.end( ), '\n' ), string.end( ) );
}

//! Convert SGP4 ECI object to state vector.
Vector6 getStateVector( const Eci state )
{
    Vector6 result( 6 );
    result[ 0 ] = state.Position( ).x;
    result[ 1 ] = state.Position( ).y;
    result[ 2 ] = state.Position( ).z;
    result[ 3 ] = state.Velocity( ).x;
    result[ 4 ] = state.Velocity( ).y;
    result[ 5 ] = state.Velocity( ).z;
    return result;
}

Vector6 getStateVectorInMetre( const Eci state )
{
    Vector6 result( 6 );
    result[ 0 ] = state.Position( ).x * 1000.0;
    result[ 1 ] = state.Position( ).y * 1000.0;
    result[ 2 ] = state.Position( ).z * 1000.0;
    result[ 3 ] = state.Velocity( ).x * 1000.0;
    result[ 4 ] = state.Velocity( ).y * 1000.0;
    result[ 5 ] = state.Velocity( ).z * 1000.0;
    return result;
}

Vector6 standardKeplerian( const Vector6 elements )
{
    Vector6 result( 6 );
    result[ 0 ] = elements[ 0 ] / 1000.0;
    result[ 1 ] = elements[ 1 ];
    result[ 2 ] = sml::convertRadiansToDegrees( elements[ 2 ] );
    result[ 3 ] = sml::convertRadiansToDegrees( elements[ 3 ] );
    result[ 4 ] = sml::convertRadiansToDegrees( elements[ 4 ] );
    result[ 5 ] = sml::convertRadiansToDegrees( elements[ 5 ] );
    return result;
}


int main( void )
{
    // conversion from km to m
    // const double km2m = 1000; 
    const double Rearth = kXKMPER; // earth radius in km
    // earth radius and diameter
    // const double EarthRadius = kXKMPER * km2m; // unit m
    // const double EarthDiam = 2 * EarthRadius;
    
    // grav. parameter 'mu' of earth
    const double muEarth = kMU*( pow( 10, 9 ) ); // unit m^3/s^2
    
    // vectors to store arrival and departure velocities for the transfer trajectory
    Vector3 DepartureVelocity( 3 );
    Vector3 ArrivalVelocity( 3 );

    // read the TLE file. Line based parsing, using string streams
    std::string line;
    
    // std::ifstream tlefile( "../../src/napa_prograde_catalog.txt" );
    // const bool is_retro = false;

    std::ifstream tlefile( "../../src/napa_retrograde_catalog.txt" );
    const bool is_retro = true;
    
    if( !tlefile.is_open() )
        perror("error while opening file");

    std::vector < Tle > tleObjects; // vector of TLE objects
    

    while( !tlefile.eof( ) )
    {
        std::vector < std::string > tleStrings;

        std::getline( tlefile, line ); // read line from the catalog file
        removeNewline( line ); // remove new line characters such as '/n'
        tleStrings.push_back( line );
        // std::cout << line << std::endl;

        std::getline( tlefile, line );
        removeNewline( line );
        tleStrings.push_back( line );
        // std::cout << line << std::endl;

        std::getline( tlefile, line );
        removeNewline( line );
        tleStrings.push_back( line );
        // std::cout << line << std::endl << std::endl;

        tleObjects.push_back( Tle( tleStrings[ 0 ], tleStrings[ 1 ], tleStrings[ 2 ] ) );
    }
    tlefile.close( );
    
    const int DebrisObjects = tleObjects.size( );
    std::cout << "Total debris objects = " << DebrisObjects << std::endl; 
    
    // some variables for the "catch" segment
    int failCount = 0;
    int catchDepartureID;
    int catchArrivalID;

    std::cout.precision( 15 );
    
    //*********************** Inputs *********************************************************************************************************//
    Tle departureObject = tleObjects[ 2 ]; 
    Tle arrivalObject = tleObjects[ 3 ];    
    DateTime departureEpoch;
    // departureEpoch.Initialise( 2016, 1, 11, 4, 51, 21, 455936 );
    departureEpoch.Initialise( 2016, 1, 11, 4, 50, 21, 0 );
    const double TOF = 15130.0; // time of flight                     
    //***************************************************************************************************************************************//
    
    SGP4 sgp4Departure( departureObject );
    const Eci tleDepartureState = sgp4Departure.FindPosition( departureEpoch );
    const Vector6 departureState = getStateVector( tleDepartureState );
    
    array3 departurePosition;
    array3 departureVelocity;
    for( int j = 0; j < 3; j++ )
    {
        departurePosition[ j ] = departureState[ j ];
        departureVelocity[ j ] = departureState[ j + 3 ];
    }
    
    const int departureObjectId = static_cast< int >( departureObject.NoradNumber( ) );
    catchDepartureID = departureObjectId; // for the catch segment of the program
                    
                    
    SGP4 sgp4Arrival( arrivalObject );
    const int arrivalObjectId = static_cast< int >( arrivalObject.NoradNumber( ) );
    catchArrivalID = arrivalObjectId;

    const DateTime arrivalEpoch = departureEpoch.AddSeconds( TOF );
    const Eci tleArrivalState = sgp4Arrival.FindPosition( arrivalEpoch );
    // const Eci tleArrivalState = sgp4Arrival.FindPosition( 0.0 );
    const Vector6 arrivalState = getStateVector( tleArrivalState );

    array3 arrivalPosition;
    array3 arrivalVelocity;
    for( int j = 0; j < 3; j++ )
    {
        arrivalPosition[ j ] = arrivalState[ j ];
        arrivalVelocity[ j ] = arrivalState[ j + 3 ];
    } 

    kep_toolbox::lambert_problem targeter( departurePosition, arrivalPosition, TOF, kMU, is_retro, 5 );
    const int numberOfSolutions = targeter.get_v1( ).size( );
    std::vector< array3 > departureDeltaVs( numberOfSolutions ); // delta-V components at the departure point
    std::vector< array3 > arrivalDeltaVs( numberOfSolutions );                
    std::vector< Real > transferDeltaVs( numberOfSolutions ); // magnitude of the total delta-V of one transfer between two points

    for ( int j = 0; j < numberOfSolutions; j++ )
    {
        array3 transferDepartureVelocity = targeter.get_v1( )[ j ]; // velocity of the s/c at the departure point in the transfer orbit
        array3 transferArrivalVelocity = targeter.get_v2( )[ j ];

        departureDeltaVs[ j ] = sml::add( transferDepartureVelocity, sml::multiply( departureVelocity, -1.0 ) );
        arrivalDeltaVs[ j ] = sml::add( transferArrivalVelocity, sml::multiply( arrivalVelocity, -1.0 ) );

        transferDeltaVs[ j ] = sml::norm< Real >( departureDeltaVs[ j ] ) + sml::norm< Real >( arrivalDeltaVs[ j ] );
    }

    const std::vector< Real >::iterator minDeltaVIterator = std::min_element( transferDeltaVs.begin( ), transferDeltaVs.end( ) );
    const int minimumDeltaVIndex = std::distance( transferDeltaVs.begin( ), minDeltaVIterator );

    Real lambertDepartureBurn = sml::norm< Real >( departureDeltaVs[ minimumDeltaVIndex ] );
    Real lambertArrivalBurn = sml::norm< Real >( arrivalDeltaVs[ minimumDeltaVIndex ] );
    Real lambertDV = transferDeltaVs[ minimumDeltaVIndex ];

    array3 minIndexDepartureVelocity = targeter.get_v1( )[ minimumDeltaVIndex ]; // best guess for velocity in transfer orbit at the departure point
    
    Vector6 lambertDepState( 6 );
    for ( int i = 0; i < 3; i++ )
    {
        lambertDepState[ i ] = departurePosition [ i ] * 1000.0;
        lambertDepState[ i + 3 ] = minIndexDepartureVelocity[ i ] * 1000.0;
    }
    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( );
    Vector6 LambertKep( 6 );
    LambertKep = astro::convertCartesianToKeplerianElements( lambertDepState, muEarth, tolerance );
    Vector6 LambertKepStd = standardKeplerian( LambertKep );
    std::cout << std::endl << "LambertKep = " << LambertKepStd << std::endl << std::endl;

    Vector3 departureVelocityGuess( 3 );
    Vector3 atomDeparturePosition( 3 );
    Vector3 atomArrivalPosition( 3 );

    for( int j = 0; j < 3; j++ )
    {
        departureVelocityGuess[ j ] = minIndexDepartureVelocity[ j ];
        atomDeparturePosition[ j ] = departurePosition[ j ];
        atomArrivalPosition[ j ] = arrivalPosition[ j ];
    }

    array3 atomDepartureVelocity;
    array3 atomArrivalVelocity;
                            
    std::string SolverStatusSummary;
    int numberOfIterations;
    const int maxIterations = 100;
    const Tle referenceTle = Tle( );
                            
    try
    {
        Vector3 outputDepartureVelocity( 3 );
        Vector3 outputArrivalVelocity( 3 );
        atom::executeAtomSolver< Real, Vector3 >( atomDeparturePosition, 
                                                  departureEpoch, 
                                                  atomArrivalPosition, 
                                                  TOF, 
                                                  departureVelocityGuess,
                                                  outputDepartureVelocity,
                                                  outputArrivalVelocity, 
                                                  SolverStatusSummary, 
                                                  numberOfIterations, 
                                                  referenceTle, 
                                                  kMU, 
                                                  kXKMPER, 
                                                  1.0e-10, 
                                                  1.0e-5, 
                                                  maxIterations );
                            
        
        Vector6 atomDepState( 6 );
        for ( int i = 0; i < 3; i++ )
        {
            atomDepState[ i ] = departurePosition[ i ] * 1000.0;
            atomDepState[ i + 3 ] = outputDepartureVelocity[ i ] * 1000.0;
        }
        Vector6 atomDepKep = astro::convertCartesianToKeplerianElements( atomDepState, muEarth, tolerance );
        Vector6 atomDepKepStd( 6 );
        atomDepKepStd = standardKeplerian( atomDepKep );
        std::cout << std::endl << "atom dep kep = " << atomDepKepStd << std::endl << std::endl;

        for( int k = 0; k < 3; k++)
        {
            atomDepartureVelocity[ k ] = outputDepartureVelocity[ k ];
            atomArrivalVelocity[ k ] = outputArrivalVelocity[ k ];
        }
                                
        array3 atomDepartureDeltaV;
        array3 atomArrivalDeltaV;
        Real AtomDeltaV;

        atomDepartureDeltaV = sml::add( atomDepartureVelocity, sml::multiply( departureVelocity, -1.0 ) );
        atomArrivalDeltaV = sml::add( atomArrivalVelocity, sml::multiply( arrivalVelocity, -1.0 ) );

        double AtomDepartureBurn = sml::norm< Real >( atomDepartureDeltaV );
        double AtomArrivalBurn = sml::norm< Real >( atomArrivalDeltaV );
        AtomDeltaV = sml::norm< Real >( atomDepartureDeltaV ) + sml::norm< Real >( atomArrivalDeltaV );
        std::cout << "Atom DeltaV = " << AtomDeltaV << std::endl;
    }
    
    catch( const std::exception& err )   
    {
        ++failCount;
        std::cout << "Exception Caught = " << err.what( ) << std::endl;
        std::cout << "For departure ID = " << catchDepartureID << " " << "For Arrival ID = " << catchArrivalID << std::endl << std::endl;
        // std::cout << "Fail count = " << failCount << std::endl;
    }

   return EXIT_SUCCESS;
}