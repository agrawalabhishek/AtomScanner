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
    
    std::ifstream tlefile( "../../src/napa_prograde_catalog.txt" );
    const bool is_retro = false;
    
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

    std::ofstream ephemerisfile;
    ephemerisfile.precision( 15 );
    std::ofstream osculatingfile;
    osculatingfile.precision( 15 );
    
    std::cout.precision( 15 );
    
    //*********************** Inputs *********************************************************************************************************//
    Tle departureObject = tleObjects[ 0 ]; 
    Tle arrivalObject = tleObjects[ 1 ];    
    DateTime departureEpoch = DateTime( 2016, 1, 13, 22, 06, 45 );
    const double TOF = 24910.0; // time of flight                     
    ephemerisfile.open( "/home/abhishek/Dropbox/Dinamica Internship/16-03-01 Napa osculating elements and ephemeris/16-03-01 ephemeris transfer case 5.csv" );
    osculatingfile.open( "/home/abhishek/Dropbox/Dinamica Internship/16-03-01 Napa osculating elements and ephemeris/16-03-01 osculating transfer case 5.csv" );
    //***************************************************************************************************************************************//
    

    SGP4 sgp4Departure( departureObject );
    const Eci tleDepartureState = sgp4Departure.FindPosition( departureEpoch );
    // const Eci tleDepartureState = sgp4Departure.FindPosition( 0.0 );
    const Vector6 departureState = getStateVector( tleDepartureState );
    
    array3 departurePosition;
    array3 departureVelocity;
    for( int j = 0; j < 3; j++ )
    {
        departurePosition[ j ] = departureState[ j ];
        departureVelocity[ j ] = departureState[ j + 3 ];
    }
    
    const int departureObjectId = static_cast< int >( departureObject.NoradNumber( ) );
    // std::cout << "departure Object ID = " << departureObjectId << std::endl;
    catchDepartureID = departureObjectId; // for the catch segment of the program
                    
                    
    SGP4 sgp4Arrival( arrivalObject );
    const int arrivalObjectId = static_cast< int >( arrivalObject.NoradNumber( ) );
    // std::cout << "Arrival Object ID = " << arrivalObjectId << std::endl;
    catchArrivalID = arrivalObjectId;

                        
    
    // DateTime arrivalEpoch = arrivalObject.Epoch( ); // take epoch from arrival TLE and then add TOF to it
    // arrivalEpoch = arrivalEpoch.AddSeconds( TOF );
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
    Vector3 departureVelocityGuess( 3 );
    Vector3 atomDeparturePosition( 3 );
    Vector3 atomArrivalPosition( 3 );
    
    Vector6 LambertState( 6 );
    array3 minLambertDepVel = targeter.get_v1( )[ minimumDeltaVIndex ];                    
    for( int i = 0; i < 3; i++ )
    {
        LambertState[ i ] = departurePosition[ i ] * 1000.0; // metre
        LambertState[ i + 3 ] = minLambertDepVel[ i ] * 1000.0; // metre/sec
    }

    array3 LambertPropPosition;
    array3 LambertPropVelocity;
    for ( int i = 0; i < 3; i++ )
    {
        LambertPropPosition[ i ] = LambertState[ i ] / 1000.0; // km
        LambertPropVelocity[ i ] = LambertState[ i + 3 ] / 1000.0; // km/s
    }

    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( );
    Vector6 LambertKep = astro::convertCartesianToKeplerianElements( LambertState, muEarth, tolerance );

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
                            
        
        Vector6 atomDepartureState( 6 );
        Vector6 atomArrivalState( 6 );

        for( int i = 0; i < 3; i++ )
        {
            atomDepartureState[ i ] = atomDeparturePosition[ i ];
            atomDepartureState[ i + 3 ] = outputDepartureVelocity[ i ];
            atomArrivalState[ i ] = atomArrivalPosition[ i ];
            atomArrivalState[ i + 3 ] = outputArrivalVelocity[ i ];
        }

        Tle transferTLE = atom::convertCartesianStateToTwoLineElements< Real, Vector6>( atomDepartureState, departureEpoch );

        // std::cout << "Transfer orbit TLE:" << testDepartureTLE << std::endl << std::endl;
        // std::cout << "Departure object TLE:" << departureObject << std::endl << std::endl;

        // Generate Ephemeris for the transfer orbit obtained from the ATOM solver
        Real TimeStep = TOF/1000.0;
        Real EphemerisJD = 0.0;
        SGP4 sgp4Ephemeris( transferTLE );
        Real tsince = 0.0;

        
        ephemerisfile << "jd" << "," << "x" << "," << "y" << "," << "z" << "," << "xdot" << "," << "ydot" << "," << "zdot" << ",";
        ephemerisfile << "lx" << "," << "ly" << "," << "lz" << "," << "lxdot" << "," << "lydot" << "," << "lzdot" << std::endl;
        
        osculatingfile << "jd" << "," << "a" << "," << "e" << "," << "i" << "," << "aop" << "," << "raan" << "," << "TA" << ",";
        osculatingfile << "raan_dot_moon" << "," << "raan_dot_sun" << "," << "raan_dot_3b" << "," << "raan_dot_j2" << "," << "raan_dot_total" << "," << "aop_dot_moon" << ",";
        osculatingfile << "aop_dot_sun" << "," << "aop_dot_3b" << "," << "aop_dot_j2" << "," << "aop_dot_total" << ",";
        osculatingfile << "La" << "," << "Le" << "," << "Li" << "," << "Laop" << "," << "Lraan" << "," << "LTA" << std::endl; 

        
        const Real j2 = 0.00108263;
        Vector6 atomDepartureStateMetre( 6 );

        for( int i = 0; i < 6; i++ )
            atomDepartureStateMetre[ i ] = atomDepartureState[ i ] * 1000;

        Vector6 nominalTransferKeplerian = astro::convertCartesianToKeplerianElements( atomDepartureStateMetre, muEarth, tolerance );
        double p_a = nominalTransferKeplerian[ 0 ];
        double p_a_km = p_a / 1000.0;
        std::cout << "semi major axis of transfer orbit = " << p_a / 1000.0 << std::endl;
        double p_e = nominalTransferKeplerian[ 1 ];
        double p_i = nominalTransferKeplerian[ 2 ];
        std::cout << "inclination = " << sml::convertRadiansToDegrees( p_i ) << std::endl;
        double p_T = 2 * sml::SML_PI * std::sqrt( std::pow( p_a, 3 ) / muEarth );
        std::cout << "time period for transfer orbit = " << p_T / 60 << std::endl;
        double p_n = ( 24 * 60 * 60 ) / p_T;
        std::cout << "mean motion for transfer orbit = " << p_n << std::endl;

        double p_aop = sml::convertRadiansToDegrees( nominalTransferKeplerian[ 3 ] );
        double p_raan = sml::convertRadiansToDegrees( nominalTransferKeplerian[ 4 ] );

        double raan_dot_moon = -0.00338 * ( std::cos( p_i ) / p_n ) / 86400.0;
        double raan_dot_sun = -0.00154 * ( std::cos( p_i ) / p_n ) / 86400.0;
        double raan_dot_j2 = -1.5 * p_n * 360.0 * j2 * std::pow( ( Rearth / p_a_km ), 2 ) * std::cos( p_i ) / ( std::pow( (1 - std::pow( p_e, 2 ) ), 2 ) * 86400.0 );
        double raan_dot_total = raan_dot_moon + raan_dot_sun + raan_dot_j2;
        std::cout << "raan_dot_total = " << raan_dot_total * 86400.0 << std::endl;

        double aop_dot_moon = 0.00169 * ( 4 - 5 * ( std::pow( std::sin( p_i ), 2 ) ) ) / ( p_n * 86400.0 );
        double aop_dot_sun = 0.00077 * ( 4 - 5 * ( std::pow( std::sin( p_i ), 2 ) ) ) / ( p_n * 86400.0 );
        double aop_dot_j2 = 0.75 * p_n * 360.0 * j2 * std::pow( ( Rearth / p_a_km ), 2 ) * ( 4 - 5 * ( std::pow( std::sin( p_i ), 2 ) ) ) 
                                                                        / ( std::pow( (1 - std::pow( p_e, 2 ) ), 2 ) * 86400.0 );
        double aop_dot_total = aop_dot_moon + aop_dot_sun + aop_dot_j2;
        std::cout << "aop_dot_total = " << aop_dot_total * 86400.0 << std::endl;
                                                               
        for( int i = 0; i < 1001; i++ )
        {
            Real t = i;
            tsince = t * TimeStep;
            DateTime EphemerisEpoch = departureEpoch; // departure epoch for the transfer orbit arc
            EphemerisEpoch = EphemerisEpoch.AddSeconds( tsince ); 

            Eci AtomEphemeris = sgp4Ephemeris.FindPosition( EphemerisEpoch );
                         
            EphemerisJD = EphemerisEpoch.ToJulian( ); // convert to Julian date

            Vector6 cart = getStateVectorInMetre( AtomEphemeris );
            
            //convert cartesian to keplerian elements using astro ( units are m and rads )
            Vector6 osc = astro::convertCartesianToKeplerianElements( cart, muEarth, tolerance ); 
                                                                       

            osculatingfile << EphemerisJD << ",";
            osculatingfile << osc[ 0 ] / 1000.0 << ",";
            osculatingfile << osc[ 1 ] << ",";
            osculatingfile << sml::convertRadiansToDegrees( osc[ 2 ] ) << ",";
            osculatingfile << sml::convertRadiansToDegrees( osc[ 3 ] ) << ",";
            osculatingfile << sml::convertRadiansToDegrees( osc[ 4 ] ) << ",";
            osculatingfile << sml::convertRadiansToDegrees( osc[ 5 ] ) << ",";
            osculatingfile << p_raan + raan_dot_moon * tsince << ",";
            osculatingfile << p_raan + raan_dot_sun * tsince << ",";
            osculatingfile << p_raan + ( raan_dot_sun + raan_dot_moon ) * tsince << ",";
            osculatingfile << p_raan + raan_dot_j2 * tsince << ",";
            osculatingfile << p_raan + raan_dot_total * tsince << ",";
            osculatingfile << p_aop + aop_dot_moon * tsince << ",";
            osculatingfile << p_aop + aop_dot_sun * tsince << ",";
            osculatingfile << p_aop + ( aop_dot_moon + aop_dot_sun ) * tsince << ",";
            osculatingfile << p_aop + aop_dot_j2 * tsince << ",";
            osculatingfile << p_aop + aop_dot_total * tsince << ",";
            osculatingfile << LambertKep[ 0 ] / 1000.0 << ",";
            osculatingfile << LambertKep[ 1 ] << ",";
            osculatingfile << sml::convertRadiansToDegrees( LambertKep[ 2 ] ) << ",";
            osculatingfile << sml::convertRadiansToDegrees( LambertKep[ 3 ] ) << ",";
            osculatingfile << sml::convertRadiansToDegrees( LambertKep[ 4 ] ) << ",";
            osculatingfile << sml::convertRadiansToDegrees( LambertKep[ 5 ] ) << std::endl;

            kep_toolbox::propagate_lagrangian( LambertPropPosition, LambertPropVelocity, TimeStep, kMU );
            for ( int j = 0; j < 3; j++ )
            {
                LambertState[ j ] = LambertPropPosition[ j ] * 1000.0; // metre
                LambertState[ j + 3 ] = LambertPropVelocity[ j ] * 1000.0; // metre/sec
            }
            LambertKep = astro::convertCartesianToKeplerianElements( LambertState, muEarth, tolerance );

            ephemerisfile << EphemerisJD << ",";
            ephemerisfile << AtomEphemeris.Position( ).x << ",";
            ephemerisfile << AtomEphemeris.Position( ).y << ",";
            ephemerisfile << AtomEphemeris.Position( ).z << ",";
            ephemerisfile << AtomEphemeris.Velocity( ).x << ",";
            ephemerisfile << AtomEphemeris.Velocity( ).y << ",";
            ephemerisfile << AtomEphemeris.Velocity( ).z << ",";
            ephemerisfile << LambertPropPosition[ 0 ] << ",";
            ephemerisfile << LambertPropPosition[ 1 ] << ",";
            ephemerisfile << LambertPropPosition[ 2 ] << ",";
            ephemerisfile << LambertPropVelocity[ 0 ] << ",";
            ephemerisfile << LambertPropVelocity[ 1 ] << ",";
            ephemerisfile << LambertPropVelocity[ 2 ] << std::endl;

            if( i == 1000 )
            {
                std::cout << "Difference in Actual Arrival Position and Final State taken from Trasnfer Orbit Ephemeris" << std::endl;
                std::cout << "X axis Difference = " << AtomEphemeris.Position( ).x - atomArrivalPosition[ 0 ] << std::endl;
                std::cout << "Y axis Difference = " << AtomEphemeris.Position( ).y - atomArrivalPosition[ 1 ] << std::endl;
                std::cout << "Z axis Difference = " << AtomEphemeris.Position( ).z - atomArrivalPosition[ 2 ] << std::endl << std::endl;
            }
        }

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
   
   ephemerisfile.close( );
   osculatingfile.close( ); 
   return EXIT_SUCCESS;
}