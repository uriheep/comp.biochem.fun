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

BOOST_AUTO_TEST_SUITE_END()
