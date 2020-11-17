// test suite: boostTestGeometry
//

#include <boost/test/unit_test.hpp>

#include "geometry.hpp"

BOOST_AUTO_TEST_SUITE( boostTestGeometry )

BOOST_AUTO_TEST_CASE( boostTestGeometry1 )
{
  const cbc::Geometry  geom;

  const std::size_t  numAtoms  =  geom.getNumAtoms();
  BOOST_CHECK_EQUAL( static_cast<std::size_t>( 0 ), numAtoms );

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
