// test suite: boostTestGeometry
//

#include <boost/test/unit_test.hpp>

#include "geometry.hpp"

BOOST_AUTO_TEST_SUITE( boostTestGeometry )

BOOST_AUTO_TEST_CASE( boostTestGeometry1 )
{
  const cbc::Geometry  geom;

  const unsigned&  numAtoms  =  geom.getNumAtoms();
  const unsigned&  numChains  =  geom.getNumChains();
  const unsigned&  numResidues  =  geom.getNumResidues();
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numAtoms );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numChains );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numResidues );

  const unsigned short  numAtomsInResidue  =  geom.getNumAtomsInResidue( 0 );
  const unsigned        numResiduesInChain  =  geom.getNumResiduesInChain( 0 );
  const char            residueName  =  geom.getResidueName( 0 );

  BOOST_CHECK_EQUAL( static_cast<unsigned short>( 0 ), numAtomsInResidue );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numResiduesInChain );
  BOOST_CHECK_EQUAL( static_cast<char>( '_' ), residueName );

  const unsigned short *  aPeriodicNumbers  =  geom.getPeriodicNumbers();
  double *  aCoordinates  =  geom.getCoordinates();
  double *  aCharges      =  geom.getCharges();
  const unsigned short *  aMultiplicity  =  geom.getMultiplicity();

  BOOST_CHECK_EQUAL( nullptr, aPeriodicNumbers );
  BOOST_CHECK_EQUAL( nullptr, aCoordinates );
  BOOST_CHECK_EQUAL( nullptr, aCharges );
  BOOST_CHECK_EQUAL( nullptr, aMultiplicity );
}

// test with a corrupt CIF file:
BOOST_AUTO_TEST_CASE( boostTestGeometry2 )
{
  cbc::Geometry  geom;
  geom.readCIF( "corrupt_1.cif" );
//  geom.readCIF( "hydrogen.cif" );

  const unsigned&  numAtoms  =  geom.getNumAtoms();
  const unsigned&  numChains  =  geom.getNumChains();
  const unsigned&  numResidues  =  geom.getNumResidues();
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numAtoms );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numChains );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numResidues );

  const unsigned short  numAtomsInResidue  =  geom.getNumAtomsInResidue( 0 );
  const unsigned        numResiduesInChain  =  geom.getNumResiduesInChain( 0 );
  const char            residueName  =  geom.getResidueName( 0 );

  BOOST_CHECK_EQUAL( static_cast<unsigned short>( 0 ), numAtomsInResidue );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 0 ), numResiduesInChain );
  BOOST_CHECK_EQUAL( static_cast<char>( '_' ), residueName );

  const unsigned short *  aPeriodicNumbers  =  geom.getPeriodicNumbers();
  double *  aCoordinates  =  geom.getCoordinates();
  double *  aCharges      =  geom.getCharges();
  const unsigned short *  aMultiplicity  =  geom.getMultiplicity();

  BOOST_CHECK_EQUAL( nullptr, aPeriodicNumbers );
  BOOST_CHECK_EQUAL( nullptr, aCoordinates );
  BOOST_CHECK_EQUAL( nullptr, aCharges );
  BOOST_CHECK_EQUAL( nullptr, aMultiplicity );
}


BOOST_AUTO_TEST_CASE( boostTestGeometry3 )
{
  cbc::Geometry  geom;
  geom.readCIF( "test.cif" );

  const unsigned&  numAtoms  =  geom.getNumAtoms();
  const unsigned&  numChains  =  geom.getNumChains();
  const unsigned&  numResidues  =  geom.getNumResidues();
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 264 ), numAtoms );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 2 ), numChains );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 22 ), numResidues );

  for ( unsigned iChain = 0; iChain < numChains; ++iChain )
  {
    const unsigned  numResiduesInChain  =  geom.getNumResiduesInChain( iChain );
    if ( 0 == iChain )
      BOOST_CHECK_EQUAL( static_cast<unsigned>( 18 ), numResiduesInChain );
    if ( 1 == iChain )
      BOOST_CHECK_EQUAL( static_cast<unsigned>( 4 ), numResiduesInChain );
  }

  for ( unsigned iResidue = 0; iResidue < numResidues; ++iResidue )
  {
    const unsigned short  numAtomsInResidue  =  geom.getNumAtomsInResidue( iResidue );
    const char            residueName  =  geom.getResidueName( iResidue );
    if ( 0 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 13 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'A', residueName );
    }
    if ( 1 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 18 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'N', residueName );
    }
    if ( 2 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 12 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'C', residueName );
    }
    if ( 3 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 16 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'E', residueName );
    }
    if ( 4 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 20 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'Q', residueName );
    }
    if ( 5 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 10 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'G', residueName );
    }
    if ( 6 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 10 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'H', residueName );
    }
    if ( 7 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 8 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'I', residueName );
    }
    if ( 8 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 8 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'L', residueName );
    }
    if ( 9 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 9 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'K', residueName );
    }
    if ( 10 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 8 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'M', residueName );
    }
    if ( 11 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 11 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'F', residueName );
    }
    if ( 12 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 7 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'P', residueName );
    }
    if ( 13 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 6 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'S', residueName );
    }
    if ( 14 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 6 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'T', residueName );
    }
    if ( 15 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 14 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'W', residueName );
    }
    if ( 16 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 12 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'Y', residueName );
    }
    if ( 17 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 7 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'V', residueName );
    }
    if ( 18 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 27 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'R', residueName );
    }
    if ( 19 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 8 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'D', residueName );
    }
    if ( 20 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 18 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'B', residueName );
    }
    if ( 21 == iResidue )
    {
      BOOST_CHECK_EQUAL( static_cast<unsigned short>( 16 ), numAtomsInResidue );
      BOOST_CHECK_EQUAL( 'Z', residueName );
    }
  }

  constexpr unsigned short  aPeriodicNumbersMustBe[ 264 ]  =  { 7, 6, 6, 8, 6, 1, 1, 1, 1, 1, 1, 8, 1, 7, 6, 6, 8, 6, 6, 8, 7, 1, 1, 1, 1, 1, 8, 1, 7, 1, 1, 7, 6, 6, 8, 6, 16, 1, 1, 1, 1, 1, 1, 7, 6, 6, 8, 6, 6, 6, 8, 8, 1, 1, 1, 1, 1, 1, 1, 7, 6, 6, 8, 6, 6, 6, 8, 7, 1, 1, 1, 1, 1, 8, 1, 1, 1, 1, 1, 7, 6, 6, 8, 1, 1, 1, 1, 8, 1, 7, 6, 6, 8, 6, 6, 7, 6, 6, 7, 7, 6, 6, 8, 6, 6, 6, 6, 7, 6, 6, 8, 6, 6, 6, 6, 7, 6, 6, 8, 6, 6, 6, 6, 7, 7, 6, 6, 8, 6, 6, 16, 6, 7, 6, 6, 8, 6, 6, 6, 6, 6, 6, 6, 7, 6, 6, 8, 6, 6, 6, 7, 6, 6, 8, 6, 8, 6, 6, 8, 6, 8, 6, 7, 6, 6, 8, 6, 6, 6, 6, 7, 6, 6, 6, 6, 6, 7, 6, 6, 8, 6, 6, 6, 6, 6, 6, 6, 8, 7, 6, 6, 8, 6, 6, 6, 7, 6, 6, 8, 6, 6, 6, 7, 6, 7, 7, 1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 6, 6, 8, 6, 6, 8, 8, 7, 6, 6, 8, 6, 6, 8, 7, 1, 1, 1, 1, 1, 8, 1, 7, 1, 1, 7, 6, 6, 8, 6, 6, 6, 8, 8, 1, 1, 1, 1, 1, 1, 1 };
  const unsigned short *  aPeriodicNumbers  =  geom.getPeriodicNumbers();
  double *  aCoordinates  =  geom.getCoordinates();

  for ( unsigned iAtom = 0; iAtom < numAtoms; ++iAtom )
    BOOST_CHECK_EQUAL( aPeriodicNumbersMustBe[ iAtom ], aPeriodicNumbers[ iAtom ] );

  for ( unsigned  i = 0; i < 3 * numAtoms; ++i )
    BOOST_CHECK_EQUAL( static_cast<double>( i + 1 ), aCoordinates[ i ] );

  double *  aCharges      =  geom.getCharges();
  const unsigned short *  aMultiplicity  =  geom.getMultiplicity();

  BOOST_CHECK_EQUAL( nullptr, aCharges );
  BOOST_CHECK_EQUAL( nullptr, aMultiplicity );
}

BOOST_AUTO_TEST_CASE( boostTestGeometry4 )
{
  cbc::Geometry  geom;
  geom.readCIF( "hydrogen.cif" );

  const unsigned&  numAtoms  =  geom.getNumAtoms();
  const unsigned&  numChains  =  geom.getNumChains();
  const unsigned&  numResidues  =  geom.getNumResidues();
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 1 ), numAtoms );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 1 ), numChains );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 1 ), numResidues );

  const unsigned short  numAtomsInResidue  =  geom.getNumAtomsInResidue( 0 );
  const unsigned        numResiduesInChain  =  geom.getNumResiduesInChain( 0 );
  const char            residueName  =  geom.getResidueName( 0 );

  BOOST_CHECK_EQUAL( static_cast<unsigned short>( 1 ), numAtomsInResidue );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 1 ), numResiduesInChain );
  BOOST_CHECK_EQUAL( static_cast<char>( 'A' ), residueName );

  const unsigned short *  aPeriodicNumbers  =  geom.getPeriodicNumbers();
  double *  aCoordinates  =  geom.getCoordinates();

  BOOST_CHECK_EQUAL( static_cast<unsigned short>( 1 ), aPeriodicNumbers[ 0 ] );
  BOOST_CHECK_CLOSE( static_cast<double>( 30.178 ), aCoordinates[ 0 ], 1 );
  BOOST_CHECK_CLOSE( static_cast<double>( 8.269 ), aCoordinates[ 1 ], 1 );
  BOOST_CHECK_CLOSE( static_cast<double>( 4.313 ), aCoordinates[ 2 ], 1 );

  double *  aCharges      =  geom.getCharges();
  const unsigned short *  aMultiplicity  =  geom.getMultiplicity();

  BOOST_CHECK_EQUAL( nullptr, aCharges );
  BOOST_CHECK_EQUAL( nullptr, aMultiplicity );
}

BOOST_AUTO_TEST_CASE( boostTestGeometry5 )
{
  cbc::Geometry  geom;
  geom.readCIF( "test2.cif" );

  const unsigned&  numAtoms  =  geom.getNumAtoms();
  const unsigned&  numChains  =  geom.getNumChains();
  const unsigned&  numResidues  =  geom.getNumResidues();
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 5 ), numAtoms );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 5 ), numChains );
  BOOST_CHECK_EQUAL( static_cast<unsigned>( 5 ), numResidues );

  for ( unsigned iResidue = 0; iResidue < numResidues; ++iResidue )
  {
    const unsigned short  numAtomsInResidue  =  geom.getNumAtomsInResidue( iResidue );
    const char            residueName  =  geom.getResidueName( 0 );
    BOOST_CHECK_EQUAL( static_cast<unsigned short>( 1 ), numAtomsInResidue );
    BOOST_CHECK_EQUAL( static_cast<char>( 'A' ), residueName );
  }
  for ( unsigned iChain = 0; iChain < numChains; ++iChain )
  {
    const unsigned  numResiduesInChain  =  geom.getNumResiduesInChain( iChain );
    BOOST_CHECK_EQUAL( static_cast<unsigned>( 1 ), numResiduesInChain );
  }


  const unsigned short *  aPeriodicNumbers  =  geom.getPeriodicNumbers();
  double *  aCoordinates  =  geom.getCoordinates();

  for ( unsigned iAtom = 0; iAtom < numAtoms; ++iAtom )
    BOOST_CHECK_EQUAL( static_cast<unsigned short>( 1 ), aPeriodicNumbers[ iAtom ] );

  for ( unsigned iAtom = 0; iAtom < numAtoms; iAtom += 3 )
  {
    BOOST_CHECK_CLOSE( static_cast<double>( 30.178 ), aCoordinates[ iAtom ], 1 );
    BOOST_CHECK_CLOSE( static_cast<double>( 8.269 ), aCoordinates[ iAtom + 1 ], 1 );
    BOOST_CHECK_CLOSE( static_cast<double>( 4.313 ), aCoordinates[ iAtom + 2 ], 1 );
  }

  double *  aCharges      =  geom.getCharges();
  const unsigned short *  aMultiplicity  =  geom.getMultiplicity();

  BOOST_CHECK_EQUAL( nullptr, aCharges );
  BOOST_CHECK_EQUAL( nullptr, aMultiplicity );
}



BOOST_AUTO_TEST_SUITE_END()
