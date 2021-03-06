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
#include <math.h>
#include <tgmath.h>

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

Vector3 crossProduct( Vector3 firstVec, Vector3 secondVec )
{
    Vector3 result( 3 );
    result[ 0 ] = firstVec[ 1 ] * secondVec[ 2 ] - secondVec[ 1 ] * firstVec[ 2 ];
    result[ 1 ] = -1.0 * ( firstVec[ 0 ] * secondVec[ 2 ] - secondVec[ 0 ] * firstVec[ 2 ] );
    result[ 2 ] = firstVec[ 0 ] * secondVec[ 1 ] - secondVec[ 0 ] * firstVec[ 1 ];
    return result;
}

Vector3 ECI2RTN( array3 ECI, Vector2D transform )
{
    Vector3 result( 3 );
    for ( int i = 0; i < 3; i++ )
    {
        Real sum = 0;
        for ( int j = 0; j < 3; j++ )
        {
            sum += transform[ i ][ j ] * ECI[ j ];
        }
        result[ i ] = sum;
    }
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
    
    std::ofstream outputfile;
    outputfile.precision( 15 );
    outputfile.open( "/home/abhishek/Dropbox/Dinamica Internship/16-03-03 Napa transfer summary DV-RTN.csv" );

    outputfile << "Departure ID" << "," << "Arrival ID" << "," << "Departure Epoch" << "," << "TOF [s]" << ",";
    outputfile << "Atom Departure Burn [m/s]" << "," << "Atom Arrival Burn [m/s]" << "," << "Atom total DV [m/s]" << "," << "Number of revolutions" << ",";
    outputfile << "Atom Iterations" << "," << "Lambert Departure Burn [m/s]" << "," << "Lambert Arrival Burn [m/s]" << ",";
    outputfile << "Lambert total DV [m/s]" << "," << "Atom Dep. Burn - Lambert Dep. Burn [m/s]" << ",";
    outputfile << "Atom Arr. Burn - Lambert Arr. Burn [m/s]" << "," << "Atom Total Burn - Lambert Total Burn [m/s]" << ",";
    outputfile << "Atom Vx Dep [m/s]" << "," << "Atom Vy Dep [m/s]" << "," << "Atom Vz Dep [m/s]" << ",";
    outputfile << "lambert Vx Dep [m/s]" << "," << "lambert Vy Dep [m/s]" << "," << "lambert Vz Dep [m/s]" << ",";
    outputfile << "Dep Burn Angle [deg]" << ",";
    outputfile << "Atom Vx Arr [m/s]" << "," << "Atom Vy Arr [m/s]" << "," << "Atom Vz Arr [m/s]" << ",";
    outputfile << "lambert Vx Arr [m/s]" << "," << "lambert Vy Arr [m/s]" << "," << "lambert Vz Arr [m/s]" << ",";
    outputfile << "Arr Burn Angle [deg]" << std::endl;
  
 for ( int runcase = 1; runcase < 6; runcase++ )
 {   
    //*********************** ****** *********************************************************************************************************//
    int transferCase = runcase;

    //****************************************************************************************************************************************//
    std::ifstream tlefile;
    bool is_retro;
    if ( transferCase == 5 )
    {
        tlefile.open( "../../src/napa_prograde_catalog.txt" );
        is_retro = false;
    }
    else 
    {
        tlefile.open( "../../src/napa_retrograde_catalog.txt" );
        is_retro = true;
    }

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
    
    
    //**************************************** Inputs ******************************************************************//
    Tle departureObject;
    Tle arrivalObject;
    DateTime departureEpoch;
    double TOF; // time of flight

    if ( transferCase == 1 )
    {
        departureObject = tleObjects[ 3 ]; 
        arrivalObject = tleObjects[ 2 ];    
        departureEpoch.Initialise( 2016, 1, 11, 4, 53, 22, 0 );
        TOF = 2830.0; // time of flight            
    }
    else if ( transferCase == 2 )
    {
        departureObject = tleObjects[ 3 ]; 
        arrivalObject = tleObjects[ 2 ];    
        departureEpoch.Initialise( 2016, 1, 11, 5, 50, 2, 0 );
        TOF = 1810.0; // time of flight            
    }
    else if ( transferCase == 3 )
    {
        departureObject = tleObjects[ 2 ]; 
        arrivalObject = tleObjects[ 3 ];    
        departureEpoch.Initialise( 2016, 1, 11, 4, 50, 21, 0 );
        TOF = 5530.0; // time of flight            
    }
    else if ( transferCase == 4 )
    {
        departureObject = tleObjects[ 2 ]; 
        arrivalObject = tleObjects[ 3 ];    
        departureEpoch.Initialise( 2016, 1, 11, 4, 50, 21, 0 );
        TOF = 15130.0; // time of flight                     
    }
    else if ( transferCase == 5 )
    {
        departureObject = tleObjects[ 0 ]; 
        arrivalObject = tleObjects[ 1 ];    
        departureEpoch.Initialise( 2016, 1, 13, 22, 6, 45, 0 );
        TOF = 24910.0; // time of flight                     
    }
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
    array3 minIndexDepartureVelocity = targeter.get_v1( )[ minimumDeltaVIndex ]; // best guess for velocity in transfer orbit at the departure point
    //***********************************************************************************************//
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
    std::cout << std::endl << transferCase << "." << '\t' << "LambertKep = " << LambertKepStd << std::endl << std::endl;
    //**************************************************************************************************//
    
    Real lambertDepartureBurn = sml::norm< Real >( departureDeltaVs[ minimumDeltaVIndex ] );
    Real lambertArrivalBurn = sml::norm< Real >( arrivalDeltaVs[ minimumDeltaVIndex ] );
    Real lambertDV = transferDeltaVs[ minimumDeltaVIndex ];

    array3 lambertDepartureDeltaV;
    array3 lambertArrivalDeltaV;
    lambertDepartureDeltaV = departureDeltaVs[ minimumDeltaVIndex ];
    lambertArrivalDeltaV = arrivalDeltaVs[ minimumDeltaVIndex ];
    
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
                            
        
        for( int k = 0; k < 3; k++)
        {
            atomDepartureVelocity[ k ] = outputDepartureVelocity[ k ];
            atomArrivalVelocity[ k ] = outputArrivalVelocity[ k ];
        }
    
        //************************************************************************************************//                        
        Vector6 atomDepState( 6 );
        for ( int i = 0; i < 3; i++ )
        {
            atomDepState[ i ] = departurePosition[ i ] * 1000.0;
            atomDepState[ i + 3 ] = outputDepartureVelocity[ i ] * 1000.0;
        }
        Vector6 atomDepKep = astro::convertCartesianToKeplerianElements( atomDepState, muEarth, tolerance );
        Vector6 atomDepKepStd( 6 );
        atomDepKepStd = standardKeplerian( atomDepKep );
        std::cout << std::endl << transferCase << "." << '\t' << "atom dep kep = " << atomDepKepStd << std::endl << std::endl;
        //************************************************************************************************//
        double timeperiod = 2 * sml::SML_PI * std::sqrt( std::pow( atomDepKepStd[ 0 ], 3 ) / kMU ); // seconds
        int revolutions = std::floor ( TOF / timeperiod );
        std::cout << "Number of revolutions = " << revolutions << std::endl;

        array3 atomDepartureDeltaV;
        array3 atomArrivalDeltaV;
        Real AtomDeltaV;

        atomDepartureDeltaV = sml::add( atomDepartureVelocity, sml::multiply( departureVelocity, -1.0 ) );
        atomArrivalDeltaV = sml::add( atomArrivalVelocity, sml::multiply( arrivalVelocity, -1.0 ) );

        double AtomDepartureBurn = sml::norm< Real >( atomDepartureDeltaV );
        double AtomArrivalBurn = sml::norm< Real >( atomArrivalDeltaV );
        AtomDeltaV = sml::norm< Real >( atomDepartureDeltaV ) + sml::norm< Real >( atomArrivalDeltaV );
        std::cout << "Atom DeltaV = " << AtomDeltaV << std::endl;

        //***********************ECI2RTN rotation matrix*******************************************************//
        Vector3 unitR( 3 );
        unitR[ 0 ] = departurePosition[ 0 ] / sml::norm< Real >( departurePosition ); 
        unitR[ 1 ] = departurePosition[ 1 ] / sml::norm< Real >( departurePosition ); 
        unitR[ 2 ] = departurePosition[ 2 ] / sml::norm< Real >( departurePosition ); 

        Vector3 RCrossV( 3 );
        RCrossV = crossProduct( atomDeparturePosition, outputDepartureVelocity );
        Real normRCrossV = sml::norm< Real >( RCrossV );
        Vector3 unitN( 3 );
        unitN[ 0 ] = RCrossV[ 0 ] / normRCrossV;
        unitN[ 1 ] = RCrossV[ 1 ] / normRCrossV;
        unitN[ 2 ] = RCrossV[ 2 ] / normRCrossV;

        Vector3 unitT( 3 );
        unitT = crossProduct( unitN, unitR );

        Vector2D transform( 3, std::vector< Real >( 3 ) );
        for ( int i = 0; i < 3; i++ )
        {
            transform[ 0 ][ i ] = unitR[ i ];
            transform[ 1 ][ i ] = unitT[ i ];
            transform[ 2 ][ i ] = unitN[ i ]; 
        }

        std::cout << "transform: " << std::endl;
        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                std::cout << transform[ i ][ j ] << '\t';
            }
            std::cout << std::endl;
        }
        

        Vector3 atomDepRTN( 3 );
        atomDepRTN = ECI2RTN( atomDepartureDeltaV, transform );
        // std::cout << "atomDepRTN[ x ] = " << atomDepRTN[ 0 ] << std::endl;
        // std::cout << "atomDepRTN[ y ] = " << atomDepRTN[ 1 ] << std::endl;
        // std::cout << "atomDepRTN[ z ] = " << atomDepRTN[ 2 ] << std::endl;

        Vector3 atomArrRTN( 3 );
        atomArrRTN = ECI2RTN( atomArrivalDeltaV, transform );

        Vector3 lambertDepRTN( 3 );
        lambertDepRTN = ECI2RTN( lambertDepartureDeltaV, transform );

        Vector3 lambertArrRTN( 3 );
        lambertArrRTN = ECI2RTN( lambertArrivalDeltaV, transform );
        //*************************************************************************************//

        outputfile << departureObjectId << ",";
        outputfile << arrivalObjectId << ",";
        outputfile << departureEpoch << ",";
        outputfile << TOF << ",";
        outputfile << AtomDepartureBurn * 1000.0 << ",";
        outputfile << AtomArrivalBurn * 1000.0 << ",";
        outputfile << AtomDeltaV * 1000.0 << ",";
        outputfile << revolutions << ",";
        outputfile << numberOfIterations << ",";
        outputfile << lambertDepartureBurn * 1000.0 << ",";
        outputfile << lambertArrivalBurn * 1000.0 << ",";
        outputfile << lambertDV * 1000.0 << ",";
        outputfile << (AtomDepartureBurn - lambertDepartureBurn) * 1000.0 << ",";
        outputfile << (AtomArrivalBurn - lambertArrivalBurn) * 1000.0 << ",";
        outputfile << (AtomDeltaV - lambertDV) * 1000.0 << ",";
        
        // outputfile << atomDepartureDeltaV[ 0 ] * 1000.0 << ",";
        // outputfile << atomDepartureDeltaV[ 1 ] * 1000.0 << ",";
        // outputfile << atomDepartureDeltaV[ 2 ] * 1000.0 << ",";
        // outputfile << lambertDepartureDeltaV[ 0 ] * 1000.0 << ",";
        // outputfile << lambertDepartureDeltaV[ 1 ] * 1000.0 << ",";
        // outputfile << lambertDepartureDeltaV[ 2 ] * 1000.0 << ",";
        // double Numerator = atomDepartureDeltaV[ 0 ] * lambertDepartureDeltaV[ 0 ] + atomDepartureDeltaV[ 1 ] * lambertDepartureDeltaV[ 1 ] + atomDepartureDeltaV[ 2 ] * lambertDepartureDeltaV[ 2 ];
        // double Denominator = sml::norm< Real >( atomDepartureDeltaV ) * sml::norm< Real >( lambertDepartureDeltaV );
        // double departureAngle = std::acos( Numerator / Denominator );
        // outputfile << sml::convertRadiansToDegrees( departureAngle ) << ",";

        outputfile << atomDepRTN[ 0 ] * 1000.0 << ",";
        outputfile << atomDepRTN[ 1 ] * 1000.0 << ",";
        outputfile << atomDepRTN[ 2 ] * 1000.0 << ",";
        outputfile << lambertDepRTN[ 0 ] * 1000.0 << ",";
        outputfile << lambertDepRTN[ 1 ] * 1000.0 << ",";
        outputfile << lambertDepRTN[ 2 ] * 1000.0 << ",";
        double Numerator = atomDepRTN[ 0 ] * lambertDepRTN[ 0 ] + atomDepRTN[ 1 ] * lambertDepRTN[ 1 ] + atomDepRTN[ 2 ] * lambertDepRTN[ 2 ];
        double Denominator = sml::norm< Real >( atomDepRTN ) * sml::norm< Real >( lambertDepRTN );
        double departureAngle = std::acos( Numerator / Denominator );
        outputfile << sml::convertRadiansToDegrees( departureAngle ) << ",";
        // outputfile << ",";

        outputfile << atomArrRTN[ 0 ] * 1000.0 << ",";
        outputfile << atomArrRTN[ 1 ] * 1000.0 << ",";
        outputfile << atomArrRTN[ 2 ] * 1000.0 << ",";
        outputfile << lambertArrRTN[ 0 ] * 1000.0 << ",";
        outputfile << lambertArrRTN[ 1 ] * 1000.0 << ",";
        outputfile << lambertArrRTN[ 2 ] * 1000.0 << ",";
        double Numerator2 = atomArrRTN[ 0 ] * lambertArrRTN[ 0 ] + atomArrRTN[ 1 ] * lambertArrRTN[ 1 ] + atomArrRTN[ 2 ] * lambertArrRTN[ 2 ];
        double Denominator2 = sml::norm< Real >( atomArrRTN ) * sml::norm< Real >( lambertArrRTN );
        double arrivalAngle = std::acos( Numerator2 / Denominator2 );
        outputfile << sml::convertRadiansToDegrees( arrivalAngle ) << std::endl;
        // outputfile << std::endl;
    }
    
    catch( const std::exception& err )   
    {
        ++failCount;
        std::cout << "Exception Caught = " << err.what( ) << std::endl;
        std::cout << "For departure ID = " << catchDepartureID << " " << "For Arrival ID = " << catchArrivalID << std::endl << std::endl;
        // std::cout << "Fail count = " << failCount << std::endl;
    }
 }  
   outputfile.close( );
   return EXIT_SUCCESS;
}