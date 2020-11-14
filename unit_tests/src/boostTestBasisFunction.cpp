// test suite:  testBasisFunction

#include <boost/test/unit_test.hpp>

#include "basis_set_sto-nG.hpp"
#include "grid3d.hpp"

#include <cstdio>
#include <cmath>
#include <random>


static
double
getValueHydrogen_STO3G_( const double&  xCenter,
                         const double&  yCenter,
                         const double&  zCenter,
                         const double&  x,
                         const double&  y,
                         const double&  z
                       ) noexcept;

static
double
getValueHydrogen_STO3G_HehreStewartPople_( const double&  xCenter,
                                           const double&  yCenter,
                                           const double&  zCenter,
                                           const double&  x,
                                           const double&  y,
                                           const double&  z
                                         ) noexcept;

static
double
getValueXe_STO3G_1S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_2S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_2P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_3S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_3P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_4S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_4P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_4D_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_5S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueXe_STO3G_5P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

// fictional orbital name !!! this is just to reflect
// the input file format of the JSON file.
// the way the basis functions are accounted in calculateions
// make this nomenclature insiginficant;
static
double
getValueXe_STO3G_5D_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept;

static
double
getValueC_STO3G_1S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueC_STO3G_2S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueC_STO3G_2P_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueN_STO3G_1S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueN_STO3G_2S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueN_STO3G_2P_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueO_STO3G_1S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueO_STO3G_2S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;

static
double
getValueO_STO3G_2P_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept;



BOOST_AUTO_TEST_SUITE( testBasisFunction )

BOOST_AUTO_TEST_CASE( testBasisFunction0 )
{
  const cbc::BasisSetGTO  bSet( "sto-3g.test1.json" );

  const std::size_t  numChemElements  =  bSet.getNumChemElements();
  const std::size_t  numOrbitals1     =  bSet.getNumOrbitals( 0 );
  const std::size_t  numOrbitals2     =  bSet.getNumOrbitals( 1 );
  BOOST_CHECK_EQUAL( std::size_t( 2 ), numChemElements );
  BOOST_CHECK_EQUAL( std::size_t( 1 ), numOrbitals1 );
  BOOST_CHECK_EQUAL( std::size_t( 3 ), numOrbitals2 );
}

BOOST_AUTO_TEST_CASE( testBasisFunction1 )
{
  const cbc::BasisSetGTO  bSet( "sto-3g.test1.json" );

  const std::size_t  numChemElements  =  bSet.getNumChemElements();
  const std::size_t  numOrbitals1     =  bSet.getNumOrbitals( 0 );
  const std::size_t  numOrbitals2     =  bSet.getNumOrbitals( 1 );
  BOOST_CHECK_EQUAL( std::size_t( 2 ), numChemElements );
  BOOST_CHECK_EQUAL( std::size_t( 1 ), numOrbitals1 );
  BOOST_CHECK_EQUAL( std::size_t( 3 ), numOrbitals2 );

  constexpr short   indexElement  =  0;
  constexpr short   indexOrbital  =  0;
  constexpr double  xCenter  =  0;
  constexpr double  yCenter  =  0;
  constexpr double  zCenter  =  0;
  constexpr double  x  =  0.01;
  constexpr double  y  =  0.01;
  constexpr double  z  =  0.01;

  const double  value1  =  bSet.getValue( indexElement,
                                          indexOrbital,
                                          xCenter,
                                          yCenter,
                                          zCenter,
                                          x,
                                          y,
                                          z
                                        );
  const double  value2  =  getValueHydrogen_STO3G_( xCenter,
                                                    yCenter,
                                                    zCenter,
                                                    x,
                                                    y,
                                                    z
                                                  );
  const double  value3  =  getValueHydrogen_STO3G_HehreStewartPople_( xCenter,
                                                                      yCenter,
                                                                      zCenter,
                                                                      x,
                                                                      y,
                                                                      z
                                                                    );
  printf( "val1 = %f\tval2 = %f\tval3 = %f\n", value1, value2, value3 );
}

BOOST_AUTO_TEST_CASE( testBasisFunction2 )
{
  const cbc::BasisSetGTO  bSet( "sto-3g.1.json" );

  const std::size_t  numChemElements  =  bSet.getNumChemElements();
  const std::size_t  numOrbitals     =  bSet.getNumOrbitals( 0 );
  BOOST_CHECK_EQUAL( std::size_t( 1 ), numChemElements );
  BOOST_CHECK_EQUAL( std::size_t( 11 ), numOrbitals );

//bSet.print_();

  std::random_device  rd;
  std::mt19937        gen( rd() );
  // distribution for individual spatial coordinates:
  std::uniform_real_distribution<>  dist( -50, +50 );

  constexpr std::size_t  numPointsToCheck  =  100000;
  constexpr short   indexElement  =  0;

  for ( std::size_t  iPoint = 0; iPoint < numPointsToCheck; ++iPoint )
  {
    for ( short  iOrbital = 0; (std::size_t)iOrbital < numOrbitals; ++iOrbital )
    {
      const double  xCenter  =  dist( gen );
      const double  yCenter  =  dist( gen );
      const double  zCenter  =  dist( gen );
      const double  x  =  dist( gen );
      const double  y  =  dist( gen );
      const double  z  =  dist( gen );

      const double  value1  =  bSet.getValue( indexElement,
                                              iOrbital,
                                              xCenter,
                                              yCenter,
                                              zCenter,
                                              x,
                                              y,
                                              z
                                            );
      double  value2  =  0;
      if ( 0 == iOrbital )
        value2  =  getValueXe_STO3G_1S_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 1 == iOrbital )
        value2  =  getValueXe_STO3G_2S_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 2 == iOrbital )
        value2  =  getValueXe_STO3G_2P_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 3 == iOrbital )
        value2  =  getValueXe_STO3G_3S_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 4 == iOrbital )
        value2  =  getValueXe_STO3G_3P_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 5 == iOrbital )
        value2  =  getValueXe_STO3G_4S_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 6 == iOrbital )
        value2  =  getValueXe_STO3G_4P_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 7 == iOrbital )
        value2  =  getValueXe_STO3G_4D_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 8 == iOrbital )
        value2  =  getValueXe_STO3G_5S_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 9 == iOrbital )
        value2  =  getValueXe_STO3G_5P_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      if ( 10 == iOrbital )
        value2  =  getValueXe_STO3G_5D_( xCenter,
                                         yCenter,
                                         zCenter,
                                         x,
                                         y,
                                         z
                                       );
      constexpr double  tolerance  =  1e-05;
      const bool  close  =  tolerance >= std::abs( value1 - value2 ) ? true : false;
      BOOST_CHECK( close );
      if ( false == close )
        printf( "iOrbital = %d\tvalue1 = %.5e\tvalue2 = %.5e\n", iOrbital, value1, value2 );
    } // for ( iOrbital )
  } // for ( iPoint )
//  printf( "val1 = %f\tval2 = %f\tval3 = %f\n", value1, value2, value3 );
}


BOOST_AUTO_TEST_CASE( testBasisFunction3 )
{
  const cbc::BasisSetGTO  bSet( "sto-3g.H_C_N_O.json" );

  const std::size_t  numChemElements  =  bSet.getNumChemElements();
  BOOST_CHECK_EQUAL( std::size_t( 4 ), numChemElements );

  for ( std::size_t  i = 0; i < numChemElements; ++i )
  {
    const std::size_t  numOrbitals  =  bSet.getNumOrbitals( i );
    if ( 0 == i )
      BOOST_CHECK_EQUAL( std::size_t( 1 ), numOrbitals );
    else
      BOOST_CHECK_EQUAL( std::size_t( 3 ), numOrbitals );
  }
//bSet.print_();

  std::random_device  rd;
  std::mt19937        gen( rd() );
  // distribution for individual spatial coordinates:
  std::uniform_real_distribution<>  dist( -50, +50 );

  constexpr std::size_t  numPointsToCheck  =  1000;

  for ( std::size_t  iPoint = 0; iPoint < numPointsToCheck; ++iPoint )
  {
    for ( std::size_t  iElement = 0; iElement < numChemElements; ++iElement )
    {
      const std::size_t  numOrbitals  =  bSet.getNumOrbitals( iElement );
      for ( std::size_t  iOrbital = 0; iOrbital < numOrbitals; ++iOrbital )
      {
        const double  xCenter  =  dist( gen );
        const double  yCenter  =  dist( gen );
        const double  zCenter  =  dist( gen );
        const double  x  =  dist( gen );
        const double  y  =  dist( gen );
        const double  z  =  dist( gen );

        const double  value1  =  bSet.getValue( iElement,
                                                iOrbital,
                                                xCenter,
                                                yCenter,
                                                zCenter,
                                                x,
                                                y,
                                                z
                                              );
        double  value2  =  0;
        if ( 0 == iOrbital )
        {
          if ( 0 == iElement )
            value2  =  getValueHydrogen_STO3G_( xCenter,
                                                yCenter,
                                                zCenter,
                                                x,
                                                y,
                                                z
                                              );
          if ( 1 == iElement )
            value2  =  getValueC_STO3G_1S_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
          if ( 2 == iElement )
            value2  =  getValueN_STO3G_1S_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
          if ( 3 == iElement )
            value2  =  getValueO_STO3G_1S_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
        }
        if ( 1 == iOrbital )
        {
          if ( 1 == iElement )
            value2  =  getValueC_STO3G_2S_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
          if ( 2 == iElement )
            value2  =  getValueN_STO3G_2S_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
          if ( 3 == iElement )
            value2  =  getValueO_STO3G_2S_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
        }
        if ( 2 == iOrbital )
        {
          if ( 1 == iElement )
            value2  =  getValueC_STO3G_2P_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
          if ( 2 == iElement )
            value2  =  getValueN_STO3G_2P_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
          if ( 3 == iElement )
            value2  =  getValueO_STO3G_2P_( xCenter,
                                            yCenter,
                                            zCenter,
                                            x,
                                            y,
                                            z
                                          );
        }
        constexpr double  tolerance  =  1e-05;
        const bool  close  =  tolerance >= std::abs( value1 - value2 ) ? true : false;
        BOOST_CHECK( close );
        if ( false == close )
        {
          printf( "iElement = %lu\tiOrbital = %lu\tvalue1 = %.5e\tvalue2 = %.5e\n", iElement, iOrbital, value1, value2 );
          printf( "point = %f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, xCenter, yCenter, zCenter );
        }
      } // for ( iOrbital )
    } // for ( iElement )
  } // for ( iPoint )
//  printf( "val1 = %f\tval2 = %f\tval3 = %f\n", value1, value2, value3 );
}

/*
BOOST_AUTO_TEST_CASE( testBasisFunction1 )
{
  constexpr double  lowX  =  -50;
  constexpr double  highX =  +50;
  constexpr double  lowY  =  -50;
  constexpr double  highY =  +50;
  constexpr double  lowZ  =  -50;
  constexpr double  highZ =  +50;
  constexpr double  stepX =  1;
  constexpr double  stepY =  1;
  constexpr double  stepZ =  1;
  cbc::Grid3D  grid( lowX, highX,
                     lowY, highY,
                     lowZ, highZ,
                     stepX, stepY, stepZ
                   );
  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;
  grid.nodes  =  new double [ numNodes ];

  // MAP BASIS FUNCTIONS ONTO THE GRID:
  cbc::BasisSetGTO  bSet;



  delete [] grid.nodes;
  grid.nodes  =  nullptr;
}
*/

BOOST_AUTO_TEST_SUITE_END()


static
double
getValueHydrogen_STO3G_( const double&  xCenter,
                         const double&  yCenter,
                         const double&  zCenter,
                         const double&  x,
                         const double&  y,
                         const double&  z
                       ) noexcept
{
  constexpr double  exp1  =  0.3425250914E+01;
  constexpr double  exp2  =  0.6239137298E+00;
  constexpr double  exp3  =  0.1688554040E+00;
  constexpr double  coef1  =  0.1543289673E+00;
  constexpr double  coef2  =  0.5353281423E+00;
  constexpr double  coef3  =  0.4446345422E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );

  return  result;
}


static
double
getValueHydrogen_STO3G_HehreStewartPople_( const double&  xCenter,
                                           const double&  yCenter,
                                           const double&  zCenter,
                                           const double&  x,
                                           const double&  y,
                                           const double&  z
                                         ) noexcept
{
  constexpr double  exp1  =  1.09818E-01;
  constexpr double  exp2  =  4.05771E-01;
  constexpr double  exp3  =  2.22766;
  constexpr double  coef1  =  4.44635E-01;
  constexpr double  coef2  =  5.35328E-01;
  constexpr double  coef3  =  1.54329E-01;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );
//printf( "fact1  = %f\tfact2 = %f\tfact3 = %f\n", factor1, factor2, factor3 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}


static
double
getValueXe_STO3G_1S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.6264584546E+04;
  constexpr double  exp2  =  0.1141101895E+04;
  constexpr double  exp3  =  0.3088267052E+03;
  constexpr double  coef1  =  0.1543289673E+00;
  constexpr double  coef2  =  0.5353281423E+00;
  constexpr double  coef3  =  0.4446345422E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_2S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.5551398381E+03;
  constexpr double  exp2  =  0.1290025597E+03;
  constexpr double  exp3  =  0.4195563620E+02;
  constexpr double  coef1  =  -0.9996722919E-01;
  constexpr double  coef2  =  0.3995128261E+00;
  constexpr double  coef3  =  0.7001154689E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_2P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.5551398381E+03;
  constexpr double  exp2  =  0.1290025597E+03;
  constexpr double  exp3  =  0.4195563620E+02;
  constexpr double  coef1  =  0.1559162750E+00;
  constexpr double  coef2  =  0.6076837186E+00;
  constexpr double  coef3  =  0.3919573931E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_3S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.8910101433E+00;
  constexpr double  exp2  =  0.4797538759E+00;
  constexpr double  exp3  =  0.2119157236E+00;
  constexpr double  coef1  =  -0.3842642608E+00;
  constexpr double  coef2  =  -0.1972567438E+00;
  constexpr double  coef3  =  0.1375495512E+01;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_3P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.8910101433E+00;
  constexpr double  exp2  =  0.4797538759E+00;
  constexpr double  exp3  =  0.2119157236E+00;
  constexpr double  coef1  =  -0.3481691526E+00;
  constexpr double  coef2  =  0.6290323690E+00;
  constexpr double  coef3  =  0.6662832743E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}


static
double
getValueXe_STO3G_4S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.7307773504E+02;
  constexpr double  exp2  =  0.2229103845E+02;
  constexpr double  exp3  =  0.8600575622E+01;
  constexpr double  coef1  =  -0.2277635023E+00;
  constexpr double  coef2  =  0.2175436044E+00;
  constexpr double  coef3  =  0.9166769611E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_4P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.7307773504E+02;
  constexpr double  exp2  =  0.2229103845E+02;
  constexpr double  exp3  =  0.8600575622E+01;
  constexpr double  coef1  =  0.4951511155E-02;
  constexpr double  coef2  =  0.5777664691E+00;
  constexpr double  coef3  =  0.4846460366E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_4D_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.7307773504E+02;
  constexpr double  exp2  =  0.2229103845E+02;
  constexpr double  exp3  =  0.8600575622E+01;
  constexpr double  coef1  =  0.2197679508E+00;
  constexpr double  coef2  =  0.6555473627E+00;
  constexpr double  coef3  =  0.2865732590E+00;

  constexpr double  factor1  =  std::pow( 2048 * std::pow( exp1, 7 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 2048 * std::pow( exp2, 7 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 2048 * std::pow( exp3, 7 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor4  =  std::pow( 2048 * std::pow( exp1, 7 ) / 9 / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor5  =  std::pow( 2048 * std::pow( exp2, 7 ) / 9 / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor6  =  std::pow( 2048 * std::pow( exp3, 7 ) / 9 / std::pow( M_PI, 3 ), 1. / 4 );

  const double  xx  =  ( x - xCenter ) * ( x - xCenter );
  const double  yy  =  ( y - yCenter ) * ( y - yCenter );
  const double  zz  =  ( z - zCenter ) * ( z - zCenter );

  const double  xy  =  ( x - xCenter ) * ( y - yCenter );
  const double  xz  =  ( x - xCenter ) * ( z - zCenter );
  const double  yz  =  ( y - yCenter ) * ( z - zCenter );

  const double  r2  =  xx + yy + zz;

  const double  result  =  coef1 * ( factor4 * ( xy + xz + yz ) + factor1 * ( 0.5 * ( 2 * zz - xx - yy ) + std::sqrt( 3. / 4 ) * ( xx - yy ) ) ) * std::exp( -exp1 * r2 )
                         + coef2 * ( factor5 * ( xy + xz + yz ) + factor2 * ( 0.5 * ( 2 * zz - xx - yy ) + std::sqrt( 3. / 4 ) * ( xx - yy ) ) ) * std::exp( -exp2 * r2 )
                         + coef3 * ( factor6 * ( xy + xz + yz ) + factor3 * ( 0.5 * ( 2 * zz - xx - yy ) + std::sqrt( 3. / 4 ) * ( xx - yy ) ) ) * std::exp( -exp3 * r2 );

  return  result;
}


static
double
getValueXe_STO3G_5S_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.7908728280E+01;
  constexpr double  exp2  =  0.3079617799E+01;
  constexpr double  exp3  =  0.1355655337E+01;
  constexpr double  coef1  =  -0.3306100626E+00;
  constexpr double  coef2  =  0.5761095338E-01;
  constexpr double  coef3  =  0.1115578745E+01;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_5P_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.7908728280E+01;
  constexpr double  exp2  =  0.3079617799E+01;
  constexpr double  exp3  =  0.1355655337E+01;
  constexpr double  coef1  =  -0.1283927634E+00;
  constexpr double  coef2  =  0.5852047641E+00;
  constexpr double  coef3  =  0.5439442040E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueXe_STO3G_5D_( const double&  xCenter,
                      const double&  yCenter,
                      const double&  zCenter,
                      const double&  x,
                      const double&  y,
                      const double&  z
                    ) noexcept
{
  constexpr double  exp1  =  0.7908728280E+01;
  constexpr double  exp2  =  0.3079617799E+01;
  constexpr double  exp3  =  0.1355655337E+01;
  constexpr double  coef1  =  0.1250662138E+00;
  constexpr double  coef2  =  0.6686785577E+00;
  constexpr double  coef3  =  0.3052468245E+00;

  constexpr double  factor1  =  std::pow( 2048 * std::pow( exp1, 7 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 2048 * std::pow( exp2, 7 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 2048 * std::pow( exp3, 7 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor4  =  std::pow( 2048 * std::pow( exp1, 7 ) / 9 / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor5  =  std::pow( 2048 * std::pow( exp2, 7 ) / 9 / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor6  =  std::pow( 2048 * std::pow( exp3, 7 ) / 9 / std::pow( M_PI, 3 ), 1. / 4 );

  const double  xx  =  ( x - xCenter ) * ( x - xCenter );
  const double  yy  =  ( y - yCenter ) * ( y - yCenter );
  const double  zz  =  ( z - zCenter ) * ( z - zCenter );

  const double  xy  =  ( x - xCenter ) * ( y - yCenter );
  const double  xz  =  ( x - xCenter ) * ( z - zCenter );
  const double  yz  =  ( y - yCenter ) * ( z - zCenter );

  const double  r2  =  xx + yy + zz;

  const double  result  =  coef1 * ( factor4 * ( xy + xz + yz ) + factor1 * ( 0.5 * ( 2 * zz - xx - yy ) + std::sqrt( 3. / 4 ) * ( xx - yy ) ) ) * std::exp( -exp1 * r2 )
                         + coef2 * ( factor5 * ( xy + xz + yz ) + factor2 * ( 0.5 * ( 2 * zz - xx - yy ) + std::sqrt( 3. / 4 ) * ( xx - yy ) ) ) * std::exp( -exp2 * r2 )
                         + coef3 * ( factor6 * ( xy + xz + yz ) + factor3 * ( 0.5 * ( 2 * zz - xx - yy ) + std::sqrt( 3. / 4 ) * ( xx - yy ) ) ) * std::exp( -exp3 * r2 );
  return  result;
}


static
double
getValueC_STO3G_1S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.7161683735E+02;
  constexpr double  exp2  =  0.1304509632E+02;
  constexpr double  exp3  =  0.3530512160E+01;
  constexpr double  coef1  =  0.1543289673E+00;
  constexpr double  coef2  =  0.5353281423E+00;
  constexpr double  coef3  =  0.4446345422E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}


static
double
getValueC_STO3G_2S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.2941249355E+01;
  constexpr double  exp2  =  0.6834830964E+00;
  constexpr double  exp3  =  0.2222899159E+00;
  constexpr double  coef1  =  -0.9996722919E-01;
  constexpr double  coef2  =  0.3995128261E+00;
  constexpr double  coef3  =  0.7001154689E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}


static
double
getValueC_STO3G_2P_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.2941249355E+01;
  constexpr double  exp2  =  0.6834830964E+00;
  constexpr double  exp3  =  0.2222899159E+00;
  constexpr double  coef1  =  0.1559162750E+00;
  constexpr double  coef2  =  0.6076837186E+00;
  constexpr double  coef3  =  0.3919573931E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}


static
double
getValueN_STO3G_1S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.9910616896E+02;
  constexpr double  exp2  =  0.1805231239E+02;
  constexpr double  exp3  =  0.4885660238E+01;
  constexpr double  coef1  =  0.1543289673E+00;
  constexpr double  coef2  =  0.5353281423E+00;
  constexpr double  coef3  =  0.4446345422E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueN_STO3G_2S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.3780455879E+01;
  constexpr double  exp2  =  0.8784966449E+00;
  constexpr double  exp3  =  0.2857143744E+00;
  constexpr double  coef1  =  -0.9996722919E-01;
  constexpr double  coef2  =  0.3995128261E+00;
  constexpr double  coef3  =  0.7001154689E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueN_STO3G_2P_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.3780455879E+01;
  constexpr double  exp2  =  0.8784966449E+00;
  constexpr double  exp3  =  0.2857143744E+00;
  constexpr double  coef1  =  0.1559162750E+00;
  constexpr double  coef2  =  0.6076837186E+00;
  constexpr double  coef3  =  0.3919573931E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueO_STO3G_1S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.1307093214E+03;
  constexpr double  exp2  =  0.2380886605E+02;
  constexpr double  exp3  =  0.6443608313E+01;
  constexpr double  coef1  =  0.1543289673E+00;
  constexpr double  coef2  =  0.5353281423E+00;
  constexpr double  coef3  =  0.4446345422E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueO_STO3G_2S_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.5033151319E+01;
  constexpr double  exp2  =  0.1169596125E+01;
  constexpr double  exp3  =  0.3803889600E+00;
  constexpr double  coef1  =  -0.9996722919E-01;
  constexpr double  coef2  =  0.3995128261E+00;
  constexpr double  coef3  =  0.7001154689E+00;

  constexpr double  factor1  =  std::pow( 2 * exp1 / M_PI, 3. / 4 );
  constexpr double  factor2  =  std::pow( 2 * exp2 / M_PI, 3. / 4 );
  constexpr double  factor3  =  std::pow( 2 * exp3 / M_PI, 3. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * std::exp( -exp3 * r2 );
  return  result;
}

static
double
getValueO_STO3G_2P_( const double&  xCenter,
                     const double&  yCenter,
                     const double&  zCenter,
                     const double&  x,
                     const double&  y,
                     const double&  z
                   ) noexcept
{
  constexpr double  exp1  =  0.5033151319E+01;
  constexpr double  exp2  =  0.1169596125E+01;
  constexpr double  exp3  =  0.3803889600E+00;
  constexpr double  coef1  =  0.1559162750E+00;
  constexpr double  coef2  =  0.6076837186E+00;
  constexpr double  coef3  =  0.3919573931E+00;

  constexpr double  factor1  =  std::pow( 128 * std::pow( exp1, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor2  =  std::pow( 128 * std::pow( exp2, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );
  constexpr double  factor3  =  std::pow( 128 * std::pow( exp3, 5 ) / std::pow( M_PI, 3 ), 1. / 4 );

  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
  const double  result  =  coef1 * factor1 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp1 * r2 )
                         + coef2 * factor2 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp2 * r2 )
                         + coef3 * factor3 * ( ( x - xCenter ) + ( y - yCenter ) + ( z - zCenter ) ) * std::exp( -exp3 * r2 );
  return  result;
}



