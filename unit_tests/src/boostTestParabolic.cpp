// test suite: boostTestParabolic

#include <boost/test/unit_test.hpp>

#include "grid3d.hpp"

#include <cmath>

static
double
getParabolic2D( const double&  x,
                const double&  y
              )
{
  if ( 0 == y )
    return  x * x;
  if ( 0 == x )
    return  y * y;
  return  0.5 * ( x * x + y * y );
}


BOOST_AUTO_TEST_SUITE( boostTestParabolic )

BOOST_AUTO_TEST_CASE( testParabolic1 )
{
  const double  x1  =  -5;
  const double  x2  =  +5;
  const double  y1  =  -5;
  const double  y2  =  +5;
  const double  z1  =  0;
  const double  z2  =  0;
  const double  step  =  0.1;
  cbc::Grid3D  grid( x1, x2,
                     y1, y2,
                     z1, z2,
                     step, step, step
                   );
  BOOST_CHECK_EQUAL( unsigned( 101 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 101 ), grid.numNodesY );
  BOOST_CHECK_EQUAL( unsigned( 1 ), grid.numNodesZ );
  BOOST_CHECK( 0 == grid.nodes );

  grid.nodes  =  new double [ grid.numNodesX * grid.numNodesY * grid.numNodesZ ];

  double  y  =  y1;
  std::size_t  index  =  0;
  while( y2 >= y )
  {
    double  x  =  x1;
    while( x2 >= x )
    {
      const double  value  =  getParabolic2D( x, y );
      grid.nodes[ index ]  =  value;
      ++index;
      x  +=  step;
    }
    y  +=  step;
  }

  BOOST_CHECK_EQUAL( index, grid.numNodesX * grid.numNodesY * grid.numNodesZ );

  delete [] grid.nodes;
  grid.nodes  =  0;
}

BOOST_AUTO_TEST_SUITE_END()
