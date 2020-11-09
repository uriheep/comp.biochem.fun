// test suite: boostTestParabolic

#include <boost/test/unit_test.hpp>

#include "grid3d.hpp"

#include <Eigen/Sparse>
#include <cmath>
#include <iostream>


typedef Eigen::SparseMatrix<double>  SpMatrix;
typedef Eigen::Triplet<double>       EigenTriplet;


static
double
getParabolic2D( const double&  x,
                const double&  y
              ) noexcept
{
  if ( 0 == y )
    return  x * x;
  if ( 0 == x )
    return  y * y;
  return  0.5 * ( x * x + y * y );
}


static
void
insertCoefficient( const std::size_t&          indRow,
                   const long&                 iX,
                   const long&                 iY,
                   const long&                 iZ,
                   const long&                 numNodesX,
                   const long&                 numNodesY,
                   const long&                 numNodesZ,
                   const double&               w,
                   const double                wXYZ,
                   const double * const        nodes,
                   std::vector<EigenTriplet>&  coeffs,
                   Eigen::VectorXd&            b
                 ) noexcept
{
  if ( -1 > iX
    || -1 > iY
    || -1 > iZ
    || numNodesX < iX
    || numNodesY < iY
    || numNodesZ < iZ
    || 0 >= numNodesX
    || 0 >= numNodesY
    || 0 >= numNodesZ
    || 0 == nodes
     )
  {
    return; // error
  }

  if ( -1 == iX )
    b( indRow )  -=  w * nodes[ iY * numNodesX + iZ * numNodesX * numNodesY ];
  if ( numNodesX == iX )
    b( indRow )  -=  w * nodes[ numNodesX - 1 + iY * numNodesX + iZ * numNodesX * numNodesY ];
  if ( -1 == iY )
    b( indRow )  -=  w * nodes[ iX + iZ * numNodesX * numNodesY ];
  if ( numNodesY == iY )
    b( indRow )  -=  w * nodes[ iX + ( numNodesY - 1 ) * numNodesX + iZ * numNodesX * numNodesY ];
  if ( -1 == iZ )
    b( indRow )  -=  w * nodes[ iX + iY * numNodesX ];
  if ( numNodesZ == iZ )
    b( indRow )  -=  w * nodes[ iX + iY * numNodesX + ( numNodesZ - 1 ) * numNodesX * numNodesY ];

  if ( -1 != iX
    && -1 != iY
    && -1 != iZ
    && numNodesX != iX
    && numNodesY != iY
    && numNodesZ != iZ
     )
  {
    b( indRow )  *=  wXYZ;
    const std::size_t  indCol  =  iX + iY * numNodesX + iZ * numNodesX * numNodesY;
    coeffs.push_back( EigenTriplet( indRow, indCol, w ) );
  }
}

static
void
buildProblem( std::vector<EigenTriplet>&  coeffs,
              Eigen::VectorXd&            b,
              const long&                 numNodesX,
              const long&                 numNodesY,
              const long&                 numNodesZ,
              const double * const        nodes,
              const double&               stepX,
              const double&               stepY,
              const double&               stepZ
            ) noexcept
{
  if ( 0 >= numNodesX
    || 0 >= numNodesY
    || 0 >= numNodesZ
    || 0 == nodes
    || 0 >= stepX
    || 0 >= stepY
    || 0 >= stepZ
     )
  {
    return; // error
  }

  for ( long k = 0; k < numNodesZ; ++k )
    for ( long j = 0; j < numNodesY; ++j )
      for ( long i = 0; i < numNodesX; ++i )
      {
        const double  wX  =  stepY * stepY + stepZ * stepZ;
        const double  wY  =  stepX * stepX + stepZ * stepZ;
        const double  wZ  =  stepX * stepX + stepY * stepY;
        const double  wXYZ  =  stepX * stepX + stepY * stepY + stepZ * stepZ;
        const std::size_t  indRow  =  i + j * numNodesX + k * numNodesX * numNodesY;
        insertCoefficient( indRow, i - 1, j, k, numNodesX, numNodesY, numNodesZ, -wX, wXYZ, nodes, coeffs, b );
        insertCoefficient( indRow, i + 1, j, k, numNodesX, numNodesY, numNodesZ, -wX, wXYZ, nodes, coeffs, b );
        insertCoefficient( indRow, i, j - 1, k, numNodesX, numNodesY, numNodesZ, -wY, wXYZ, nodes, coeffs, b );
        insertCoefficient( indRow, i, j + 1, k, numNodesX, numNodesY, numNodesZ, -wY, wXYZ, nodes, coeffs, b );
        insertCoefficient( indRow, i, j, k - 1, numNodesX, numNodesY, numNodesZ, -wZ, wXYZ, nodes, coeffs, b );
        insertCoefficient( indRow, i, j, k + 1, numNodesX, numNodesY, numNodesZ, -wZ, wXYZ, nodes, coeffs, b );
        insertCoefficient( indRow, i, j, k, numNodesX, numNodesY, numNodesZ, wX + wY + wZ, wXYZ, nodes, coeffs, b );
      }
}



BOOST_AUTO_TEST_SUITE( boostTestParabolic )

BOOST_AUTO_TEST_CASE( testParabolic1 )
{
  const double  x1  =  -5;
  const double  x2  =  +5;
  const double  y1  =  0;
  const double  y2  =  0;
  const double  z1  =  0;
  const double  z2  =  0;
  const double  step  =  0.1;//2;//0.1;
  cbc::Grid3D  grid( x1, x2,
                     y1, y2,
                     z1, z2,
                     step, step, step
                   );
  BOOST_CHECK_EQUAL( unsigned( 101 ), grid.numNodesX );
  BOOST_CHECK_EQUAL( unsigned( 1 ), grid.numNodesY );
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
      const double  value  =  getParabolic2D( x, 0 );
      grid.nodes[ index ]  =  value;
      ++index;
      x  +=  step;
    }
    y  +=  step;
  }

  BOOST_CHECK_EQUAL( index, grid.numNodesX * grid.numNodesY * grid.numNodesZ );

  /* ********************************* */
  // finite difference scheme:
  std::vector<EigenTriplet>  coefficients;
  Eigen::VectorXd            b( grid.numNodesX * grid.numNodesY * grid.numNodesZ );
  for ( std::size_t i = 0; i < (unsigned long)( grid.numNodesX * grid.numNodesY * grid.numNodesZ ); ++i )
    b( i )  =  grid.nodes[ i ];

  buildProblem( coefficients,
                b,
                grid.numNodesX,
                grid.numNodesY,
                grid.numNodesZ,
                grid.nodes,
                step,
                step,
                step
              );

  const std::size_t  n  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;

  SpMatrix  matrix( n, n );
  matrix.setFromTriplets( coefficients.begin(), coefficients.end() );

  Eigen::SimplicialCholesky<SpMatrix>  chol( matrix );
  Eigen::VectorXd  x  =  chol.solve( b );

  std::cout << matrix << std::endl;
  std::cout << b << std::endl;
  std::cout << std::endl << x << std::endl;

  /* ********************************* */

  delete [] grid.nodes;
  grid.nodes  =  0;
}

BOOST_AUTO_TEST_CASE( testParabolic2 )
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

  /* ********************************* */
  cbc::Grid3D  gridFunction( x1, x2,
                             y1, y2,
                             z1, z2,
                             step, step, step
                           );

  gridFunction.nodes  =  new double [ gridFunction.numNodesX * gridFunction.numNodesY * gridFunction.numNodesZ ];
  y  =  y1;
  index  =  0;
  while( y2 >= y )
  {
    double  x  =  x1;
    while( x2 >= x )
    {
      const double  value  =  1;
      gridFunction.nodes[ index ]  =  value;
      ++index;
      x  +=  step;
    }
    y  +=  step;
  }
  // finite difference scheme:
  
  delete [] gridFunction.nodes;
  gridFunction.nodes  =  0;
  /* ********************************* */

  delete [] grid.nodes;
  grid.nodes  =  0;
}

BOOST_AUTO_TEST_SUITE_END()
