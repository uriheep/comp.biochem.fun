// test suite testGrid3D

#include <boost/test/unit_test.hpp>

#include "grid3d.hpp"


BOOST_AUTO_TEST_SUITE( testGrid3D )

BOOST_AUTO_TEST_CASE( testGrid3D_1 )
{
  const cbc::Grid3D  grid;
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_2 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepY  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           stepX,
                           stepY,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 11 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 11 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 11 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_3 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepY  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x2, x1,
                           y1, y2,
                           z1, z2,
                           stepX,
                           stepY,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_4 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepY  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y2, y1,
                           z1, z2,
                           stepX,
                           stepY,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_5 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepY  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z2, z1,
                           stepX,
                           stepY,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_6 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepY  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           0,
                           stepY,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_7 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepY  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           -1,
                           stepY,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_8 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           stepX,
                           0,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_9 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepZ  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           stepX,
                           -2,
                           stepZ
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_10 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepY  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           stepX,
                           stepY,
                           0
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_CASE( testGrid3D_11 )
{
  const double  x1  =  -10;
  const double  x2  =  0;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  +10;
  const double  stepX  =  1;
  const double  stepY  =  1;
  const cbc::Grid3D  grid( x1, x2,
                           y1, y2,
                           z1, z2,
                           stepX,
                           stepY,
                           -3
                         );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 0 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );
}

BOOST_AUTO_TEST_SUITE_END()
