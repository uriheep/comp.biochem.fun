#include "scf_arma_solver.hpp"

#include "grid3d.hpp"
#include "geometry.hpp"
#include "basis_set_sto-nG.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <random>
#include <cmath>


namespace  cbc {

static
double
getMarginValue_() noexcept
{
  return  5e-2; // [ A ] = margin value for the whole spatial region mapped onto a grid
}


SCFArmaSolver::SCFArmaSolver( const Geometry * const     geometry,
                              const BasisSetGTO * const  basisSet
                            ) : geometry_( geometry ),
                                basisSet_( basisSet ),
                                numMolecularOrbitals_( 0 ),
                                numBasisFunctions_( 0 ),
                                energy_( 0 ),
                                resolution_( 0 ),
                                maxNumNodesInGrid_( 0 ),
                                aGrids_( nullptr ),
                                mCoefficients_()
{
  numMolecularOrbitals_  =  getNumMolecularOrbitals_();
  numBasisFunctions_  =  getNumBasisFunctions_();
  initializeCoefficients_();
}


SCFArmaSolver::~SCFArmaSolver() noexcept
{
  delete [] aGrids_;
  aGrids_  =  nullptr;
//  delete [] aCoefficients_;
//  aCoefficients_  =  nullptr;
}

void
SCFArmaSolver::setResolution( const double&  delta ) noexcept
{
  resolution_  =  delta;
}

void
SCFArmaSolver::setMaxNumNodesInGrid( const std::size_t&  numNodes ) noexcept
{
  maxNumNodesInGrid_  =  numNodes;
}

const double&
SCFArmaSolver::getResolution() const noexcept
{
  return  resolution_;
}

std::size_t
SCFArmaSolver::getNumNodesInGrid() const noexcept
{
  if ( nullptr == aGrids_ )
    return  0;
  const unsigned  numNodesX  =  aGrids_[ 0 ].numNodesX;
  const unsigned  numNodesY  =  aGrids_[ 0 ].numNodesY;
  const unsigned  numNodesZ  =  aGrids_[ 0 ].numNodesZ;
  return  static_cast<std::size_t>( numNodesX * numNodesY * numNodesZ );
}
/*
void
SCFArmaSolver::run( const double&  tolerance ) noexcept
{
  if ( 0 == resolution_ ) // default ( unset ) value
    resolution_  =  1e-2; // [ A ] = 1e-12 [ m ]

  // get spatial limits for the whole region:
  const std::tuple<double, double>  xRange  =  getSpatialLimitsX_();
  const std::tuple<double, double>  yRange  =  getSpatialLimitsY_();
  const std::tuple<double, double>  zRange  =  getSpatialLimitsZ_();

  const arma::cx_mat  sMatrix  =  getSMatrix_( xRange, yRange, zRange );
  const arma::cx_mat  fMatrix  =  getFMatrix_( xRange, yRange, zRange );



}

void
SCFArmaSolver::run( const std::size_t&  numIterations ) noexcept
{

}
*/
const double&
SCFArmaSolver::getEnergy() const noexcept
{
  return  energy_;
}

unsigned
SCFArmaSolver::getNumMolecularOrbitals_() const noexcept
{ // ASSUMPTION: the input JSON file with the basis set 
  // contains entries sorted in the ascending
  // order of periodic numbers:
  if ( nullptr == basisSet_ )
    return  0; // error

  const std::list<unsigned>  listNumAtomsOfEachType  =  getListNumAtomsOfEachType_();
  const std::size_t  numChemElements  =  basisSet_->getNumChemElements();
  const std::size_t  numChemElementsMustBe  =  listNumAtomsOfEachType.size();
  if ( numChemElements != numChemElementsMustBe )
    return  0; // error

  unsigned  result  =  0;
  unsigned short  iChemElement  = 0;
  for ( std::list<unsigned>::const_iterator  it = listNumAtomsOfEachType.begin();
        listNumAtomsOfEachType.end() != it;
        ++it
      )
  {
    const unsigned  numAtomsOfThisType  =  *it;
    const std::size_t  numOrbitals  =  basisSet_->getNumOrbitals( iChemElement );
    result  +=  ( numOrbitals * numAtomsOfThisType );
    ++iChemElement;
  }
  return  result;
}

std::list<unsigned>
SCFArmaSolver::getListNumAtomsOfEachType_() const noexcept
{ // ASSUMPTION: the input JSON file with the basis set 
  // contains entries sorted in the ascending
  // order of periodic numbers:
  if ( nullptr == geometry_ )
    return  std::list<unsigned>();

  const unsigned&  numAtoms  =  geometry_->getNumAtoms();
  if ( 0 == numAtoms )
    return  std::list<unsigned>();

  const unsigned short * aPeriodicNumbers  =  geometry_->getPeriodicNumbers();
  std::vector<unsigned short>  vec( aPeriodicNumbers, aPeriodicNumbers + numAtoms );
  std::sort( vec.begin(), vec.end() );

  // list of numbers of atoms of each type sorted by their respective periodic numbers:
  std::list<unsigned>  result;
  unsigned  numAtomsOfGivenType  =  1;
  unsigned short  typeOfAtom  =  vec[ 0 ];

  for ( unsigned iAtom = 1; iAtom < numAtoms; ++iAtom )
  {
    const unsigned short  typeOfAtomCurrent  =  vec[ iAtom ];
    if ( typeOfAtomCurrent == typeOfAtom )
      ++numAtomsOfGivenType;
    else
      {
        result.push_back( numAtomsOfGivenType );
        typeOfAtom  =  typeOfAtomCurrent;
        numAtomsOfGivenType  =  1;
      }
  }
  return  result;
}

unsigned short
SCFArmaSolver::getNumBasisFunctions_() const noexcept
{
  if ( nullptr == basisSet_ )
    return  0; // error

  const std::size_t&  numChemElements  =  basisSet_->getNumChemElements();
  unsigned short  result  =  0;
  for ( std::size_t  iElem = 0; iElem < numChemElements; ++iElem )
  {
    const std::size_t  numOrbitals  =  basisSet_->getNumOrbitals( iElem );
    result  +=  numOrbitals;
  }
  return  result;
}

void
SCFArmaSolver::initializeCoefficients_() noexcept
{
  if ( 0 == numMolecularOrbitals_ )
    return;

  arma::cx_dmat  result( numMolecularOrbitals_, numMolecularOrbitals_, arma::fill::zeros );
  std::random_device  rd;
  std::mt19937        gen( rd() );
  std::uniform_real_distribution<>  dist( 0, 1. / numMolecularOrbitals_ );
  for ( unsigned iRow = 0; iRow < numMolecularOrbitals_; ++iRow )
    for ( unsigned iCol = 0; iCol < numMolecularOrbitals_; ++iCol )
    {
      const double  randomCoef  =  dist( gen );
      const arma::cx_double  value  =  arma::cx_double( randomCoef, 0 );
      result( iRow, iCol )  =  value;
    }
  mCoefficients_  =  result;
}

std::tuple<double, double>
SCFArmaSolver::getSpatialLimitsX_() const noexcept
{
  if ( nullptr == geometry_ )
    return  std::tuple<double, double>( 0, 0 );

  const unsigned  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  double  lowLimit  =  std::numeric_limits<double>::max();
  double  upLimit  =  std::numeric_limits<double>::min();
  for ( unsigned  iCoord = 0; iCoord < 3 * numAtoms; iCoord += 3 )
  {
    const double  x  =  aCoordinates[ iCoord ];
    if ( x > upLimit )
      upLimit  =  x;
    if ( x < lowLimit )
      lowLimit  =  x;
  }
  return  std::make_tuple( lowLimit, upLimit );
}

std::tuple<double, double>
SCFArmaSolver::getSpatialLimitsY_() const noexcept
{
  if ( nullptr == geometry_ )
    return  std::tuple<double, double>( 0, 0 );

  const unsigned  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  double  lowLimit  =  std::numeric_limits<double>::max();
  double  upLimit  =  std::numeric_limits<double>::min();
  for ( unsigned  iCoord = 1; iCoord < 3 * numAtoms; iCoord += 3 )
  {
    const double  x  =  aCoordinates[ iCoord ];
    if ( x > upLimit )
      upLimit  =  x;
    if ( x < lowLimit )
      lowLimit  =  x;
  }
  return  std::make_tuple( lowLimit, upLimit );
}

std::tuple<double, double>
SCFArmaSolver::getSpatialLimitsZ_() const noexcept
{
  if ( nullptr == geometry_ )
    return  std::tuple<double, double>( 0, 0 );

  const unsigned  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  double  lowLimit  =  std::numeric_limits<double>::max();
  double  upLimit  =  std::numeric_limits<double>::min();
  for ( unsigned  iCoord = 2; iCoord < 3 * numAtoms; iCoord += 3 )
  {
    const double  x  =  aCoordinates[ iCoord ];
    if ( x > upLimit )
      upLimit  =  x;
    if ( x < lowLimit )
      lowLimit  =  x;
  }
  return  std::make_tuple( lowLimit, upLimit );
}

arma::sp_cx_dmat
SCFArmaSolver::getSMatrix_( const std::tuple<double, double>&  xRange,
                            const std::tuple<double, double>&  yRange,
                            const std::tuple<double, double>&  zRange
                          ) const noexcept
{
  try {

    const unsigned&  numAtoms  =  geometry_->getNumAtoms();
    const double *  aCoordinates  =  geometry_->getCoordinates();
    const unsigned short * aPeriodicNumbers  =  geometry_->getPeriodicNumbers();

    std::list<std::size_t>  listIIndices;
    std::list<std::size_t>  listJIndices;
    std::list<double>  listMatrixElements;

    // to verify that the resulting sparse matrix will have the correct dimensions:
    bool  maxRowRegistered  =  false;
    bool  maxColRegistered  =  false;
    // try to account for interaction of each orbital of each atom
    // with each orbital of each other atom
    // and between orbitals within the same atom:
    std::size_t  iMatrixIndex  =  0;
    unsigned  iAtom1  =  0;
    for ( unsigned  iCoord1 = 0; iCoord1 < 3 * numAtoms; iCoord1 += 3 )
    {
      const double  xCenter1  =  aCoordinates[ iCoord1 ];
      const double  yCenter1  =  aCoordinates[ iCoord1 + 1 ];
      const double  zCenter1  =  aCoordinates[ iCoord1 + 2 ];
      const unsigned short  periodicNumber1  =  aPeriodicNumbers[ iAtom1 ];
      const std::size_t  numOrbitals1  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber1 ) );
      for ( std::size_t  iOrbital1 = 0; iOrbital1 < numOrbitals1; ++iOrbital1 )
      {
        std::size_t  jMatrixIndex  =  0;
        unsigned  iAtom2  =  0;
        for ( unsigned  iCoord2 = 0; iCoord2 < 3 * numAtoms; iCoord2 += 3 )
        {
          const double  xCenter2  =  aCoordinates[ iCoord2 ];
          const double  yCenter2  =  aCoordinates[ iCoord2 + 1 ];
          const double  zCenter2  =  aCoordinates[ iCoord2 + 2 ];
          const unsigned short  periodicNumber2  =  aPeriodicNumbers[ iAtom2 ];
          const std::size_t  numOrbitals2  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber2 ) );
          for ( std::size_t  iOrbital2 = 0; iOrbital2 < numOrbitals2; ++iOrbital2 )
          {
            const double  matrixElement  =  getSMatrixValue_( xRange, yRange, zRange,
                                                              xCenter1, yCenter1, zCenter1,
                                                              xCenter2, yCenter2, zCenter2,
                                                              periodicNumber1,
                                                              periodicNumber2,
                                                              iOrbital1,
                                                              iOrbital2
                                                            );
            const double  epsilon  =  std::numeric_limits<double>::epsilon();
            if ( epsilon < matrixElement )
            {
              listIIndices.push_back( iMatrixIndex );
              listJIndices.push_back( jMatrixIndex );
              listMatrixElements.push_back( matrixElement );
              if ( numMolecularOrbitals_ - 1 == iMatrixIndex )
                maxRowRegistered  =  true;
              if ( numMolecularOrbitals_ - 1 == jMatrixIndex )
                maxColRegistered  =  true;
            }
            ++jMatrixIndex;
          } // for ( iOrbital2 )
          ++iAtom2;
        } // for ( iCoord2 )
        ++iMatrixIndex;
      } // for ( iOrbital1 )
      ++iAtom1;
    } // for ( iCoord1 )

    if ( false == maxRowRegistered )
    {
      listIIndices.push_back( numMolecularOrbitals_ - 1 );
      listJIndices.push_back( 0 );
      listMatrixElements.push_back( 0. );
    }
    if ( false == maxColRegistered )
    {
      listIIndices.push_back( 0 );
      listJIndices.push_back( numMolecularOrbitals_ - 1 );
      listMatrixElements.push_back( 0. );
    }
    // prepare data for the Armadillo sparce matrix constructor ( it there a better solution? ):
    const std::size_t  numNonZeroElementsMatrix  =  listIIndices.size();

    arma::cx_vec  vecMatrixValues( numNonZeroElementsMatrix, arma::fill::zeros );
    arma::umat  locationsInMatrix( 2, numNonZeroElementsMatrix, arma::fill::zeros );

    std::list<std::size_t>::iterator  itIIndices  =  listIIndices.begin();
    std::list<std::size_t>::iterator  itJIndices  =  listJIndices.begin();
    std::list<double>::iterator       itMatrixElements  =  listMatrixElements.begin();
    for ( std::size_t  i = 0;
          listIIndices.end() != itIIndices
       && listJIndices.end() != itJIndices
       && listMatrixElements.end() != itMatrixElements
       && i < numNonZeroElementsMatrix
          ;
          ++i
        )
    {
      locationsInMatrix( 0, i )  =  *itIIndices;
      locationsInMatrix( 1, i )  =  *itJIndices;
      const double  matrixElement  =  *itMatrixElements;
      vecMatrixValues( i )  =  arma::cx_double( matrixElement, 0 );
      ++itIIndices;
      ++itJIndices;
      ++itMatrixElements;
    }

    const arma::sp_cx_dmat  result( locationsInMatrix, vecMatrixValues );
    return  result;
  }
  catch( const std::bad_alloc&  ba )
  {
    printf( "SCFArmaSolver: %s\n", ba.what() );
  }

  return  arma::sp_cx_dmat();
}

double
SCFArmaSolver::getSMatrixValue_( const std::tuple<double, double>&  xRange,
                                 const std::tuple<double, double>&  yRange,
                                 const std::tuple<double, double>&  zRange,
                                 const double&                      xCenter1,
                                 const double&                      yCenter1,
                                 const double&                      zCenter1,
                                 const double&                      xCenter2,
                                 const double&                      yCenter2,
                                 const double&                      zCenter2,
                                 const unsigned short&              periodicNumber1,
                                 const unsigned short&              periodicNumber2,
                                 const std::size_t&                 iOrbital1,
                                 const std::size_t&                 iOrbital2
                               ) const noexcept
{
  const double  margin  =  getMarginValue_();
  const double  lowX  =  std::get<0>( xRange ) - margin;
  const double  highX  =  std::get<1>( xRange ) + margin;
  const double  lowY  =  std::get<0>( yRange ) - margin;
  const double  highY  =  std::get<1>( yRange ) + margin;
  const double  lowZ  =  std::get<0>( zRange ) - margin;
  const double  highZ  =  std::get<1>( zRange ) + margin;

  Grid3D  grid( lowX, highX,
                lowY, highY,
                lowZ, highZ,
                resolution_,
                resolution_,
                resolution_
              );

  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;
  double  matrixElement  =  0;
  for ( std::size_t  iNode = 0; iNode < numNodes; ++iNode )
  {
    const std::size_t  iNodeZ  =  iNode / grid.numNodesX / grid.numNodesY;
    const std::size_t  iNodeY  =  ( iNode - iNodeZ * grid.numNodesX * grid.numNodesY ) / grid.numNodesX;
    const std::size_t  iNodeX  =  iNode - iNodeZ * grid.numNodesX * grid.numNodesY - iNodeY * grid.numNodesX;
    const double  x  =  lowX + iNodeX * resolution_;
    const double  y  =  lowY + iNodeY * resolution_;
    const double  z  =  lowZ + iNodeZ * resolution_;
    const double  val1  =  basisSet_->getValue( static_cast<short>( periodicNumber1 ),
                                                static_cast<short>( iOrbital1 ),
                                                xCenter1,
                                                yCenter1,
                                                zCenter1,
                                                x,
                                                y,
                                                z
                                              );
    const double  val2  =  basisSet_->getValue( static_cast<short>( periodicNumber2 ),
                                                static_cast<short>( iOrbital2 ),
                                                xCenter2,
                                                yCenter2,
                                                zCenter2,
                                                x,
                                                y,
                                                z
                                              );
    matrixElement  +=  ( val1 * val2 );
  } // for ( iNode )
  return  matrixElement;
}

arma::sp_cx_dmat
SCFArmaSolver::getLaplaceMatrix_( const std::tuple<double, double>&  xRange,
                                  const std::tuple<double, double>&  yRange,
                                  const std::tuple<double, double>&  zRange
                                ) const noexcept
{
  const double  margin  =  getMarginValue_();
  const double  lowX  =  std::get<0>( xRange ) - margin;
  const double  highX  =  std::get<1>( xRange ) + margin;
  const double  lowY  =  std::get<0>( yRange ) - margin;
  const double  highY  =  std::get<1>( yRange ) + margin;
  const double  lowZ  =  std::get<0>( zRange ) - margin;
  const double  highZ  =  std::get<1>( zRange ) + margin;

  Grid3D  grid( lowX, highX,
                lowY, highY,
                lowZ, highZ,
                resolution_,
                resolution_,
                resolution_
              );
  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;
  const unsigned  numNodesX  =  grid.numNodesX;
  const unsigned  numNodesY  =  grid.numNodesY;
  // non-zero elements of the matrix:
  // main diagonal = numNodes
  // two secondary diagonals  =  2 * ( numNodes - 2 )
  // two secondary diagonals  =  2 * ( numNodes - numNodesX - numNodesX )
  // two secondary diagonals  =  2 * ( numNodes - numNodesX * numNodesY - numNodesX * numNodesY )
  const std::size_t  numNonZeroElementsMatrix  =  numNodes + 2 * ( numNodes - 2 ) + 2 * ( numNodes - 2 * numNodesX ) + 2 * ( numNodes - 2 * numNodesX * numNodesY );

  arma::cx_vec  vecMatrixValues( numNonZeroElementsMatrix, arma::fill::zeros );
  arma::umat  locationsInMatrix( 2, numNonZeroElementsMatrix, arma::fill::zeros );

  std::size_t  matrixIndex  =  0;
  // main diagonal elements:
  for ( std::size_t  i = 0; i < numNodes; ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( -6 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex;
    locationsInMatrix( 1, i )  =  matrixIndex;
    ++matrixIndex;
  }
  // secondary upper-diagonal elements:
  matrixIndex  =  0;
  for ( std::size_t  i = numNodes; i < numNodes + ( numNodes - 2 ); ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( 1 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex;
    locationsInMatrix( 1, i )  =  matrixIndex + 1;
    ++matrixIndex;
  }
  // secondary lower-diagonal elements:
  matrixIndex  =  0;
  for ( std::size_t  i = numNodes + ( numNodes - 2 ); i < numNodes + 2 * ( numNodes - 2 ); ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( 1 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex + 1;
    locationsInMatrix( 1, i )  =  matrixIndex;
    ++matrixIndex;
  }
  // tertiary upper-diagonal elements:
  matrixIndex  =  0;
  for ( std::size_t  i = numNodes + 2 * ( numNodes - 2 ); i < numNodes + 2 * ( numNodes - 2 ) + ( numNodes - 2 * numNodesX ); ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( 1 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex;
    locationsInMatrix( 1, i )  =  matrixIndex + numNodesX;
    ++matrixIndex;
  }
  // tertiary lower-diagonal elements:
  matrixIndex  =  0;
  for ( std::size_t  i = numNodes + 2 * ( numNodes - 2 ) + ( numNodes - 2 * numNodesX ); i < numNodes + 2 * ( numNodes - 2 ) + 2 * ( numNodes - 2 * numNodesX ); ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( 1 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex + numNodesX;
    locationsInMatrix( 1, i )  =  matrixIndex;
    ++matrixIndex;
  }
  // quaternary upper-diagonal elements:
  matrixIndex  =  0;
  for ( std::size_t  i = numNodes + 2 * ( numNodes - 2 ) + 2 * ( numNodes - 2 * numNodesX ); i < numNodes + 2 * ( numNodes - 2 ) + 2 * ( numNodes - 2 * numNodesX ) + ( numNodes - 2 * numNodesX * numNodesY ); ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( 1 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex;
    locationsInMatrix( 1, i )  =  matrixIndex + numNodesX * numNodesY;
    ++matrixIndex;
  }
  // quaternary lower-diagonal elements:
  matrixIndex  =  0;
  for ( std::size_t  i = numNodes + 2 * ( numNodes - 2 ) + 2 * ( numNodes - 2 * numNodesX ) + ( numNodes - 2 * numNodesX * numNodesY ); i < numNodes + 2 * ( numNodes - 2 ) + 2 * ( numNodes - 2 * numNodesX ) + 2 * ( numNodes - 2 * numNodesX * numNodesY ); ++i )
  {
    vecMatrixValues( i )  =  arma::cx_double( 1 / resolution_ / resolution_, 0 );
    locationsInMatrix( 0, i )  =  matrixIndex + numNodesX * numNodesY;
    locationsInMatrix( 1, i )  =  matrixIndex;
    ++matrixIndex;
  }

  const arma::sp_cx_dmat  result( locationsInMatrix, vecMatrixValues );
  return  result;
}

double
SCFArmaSolver::getFMatrixKineticValue_( const arma::sp_cx_dmat&            laplaceMatrix,
                                        const std::tuple<double, double>&  xRange,
                                        const std::tuple<double, double>&  yRange,
                                        const std::tuple<double, double>&  zRange,
                                        const double&                      xCenter1,
                                        const double&                      yCenter1,
                                        const double&                      zCenter1,
                                        const double&                      xCenter2,
                                        const double&                      yCenter2,
                                        const double&                      zCenter2,
                                        const unsigned short&              periodicNumber1,
                                        const unsigned short&              periodicNumber2,
                                        const std::size_t&                 iOrbital1,
                                        const std::size_t&                 iOrbital2
                                      ) const noexcept
{
  const double  margin  =  getMarginValue_();
  const double  lowX  =  std::get<0>( xRange ) - margin;
  const double  highX  =  std::get<1>( xRange ) + margin;
  const double  lowY  =  std::get<0>( yRange ) - margin;
  const double  highY  =  std::get<1>( yRange ) + margin;
  const double  lowZ  =  std::get<0>( zRange ) - margin;
  const double  highZ  =  std::get<1>( zRange ) + margin;

  Grid3D  grid( lowX, highX,
                lowY, highY,
                lowZ, highZ,
                resolution_,
                resolution_,
                resolution_
              );

  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;

  arma::cx_vec  vecRight( numNodes, arma::fill::zeros );
  arma::cx_vec  vecLeft( numNodes, arma::fill::zeros );

  for ( std::size_t  iNode = 0; iNode < numNodes; ++iNode )
  {
    const std::size_t  iNodeZ  =  iNode / grid.numNodesX / grid.numNodesY;
    const std::size_t  iNodeY  =  ( iNode - iNodeZ * grid.numNodesX * grid.numNodesY ) / grid.numNodesX;
    const std::size_t  iNodeX  =  iNode - iNodeZ * grid.numNodesX * grid.numNodesY - iNodeY * grid.numNodesX;
    const double  x  =  lowX + iNodeX * resolution_;
    const double  y  =  lowY + iNodeY * resolution_;
    const double  z  =  lowZ + iNodeZ * resolution_;
    const double  val1  =  basisSet_->getValue( static_cast<short>( periodicNumber1 ),
                                                static_cast<short>( iOrbital1 ),
                                                xCenter1,
                                                yCenter1,
                                                zCenter1,
                                                x,
                                                y,
                                                z
                                              );
    vecRight( iNode )  =  arma::cx_double( val1, 0 );
    const double  val2  =  basisSet_->getValue( static_cast<short>( periodicNumber2 ),
                                                static_cast<short>( iOrbital2 ),
                                                xCenter2,
                                                yCenter2,
                                                zCenter2,
                                                x,
                                                y,
                                                z
                                              );
    vecLeft( iNode )  =  arma::cx_double( val2, 0 );
  } // for ( iNode )
  const arma::cx_double  resultComplex  =  arma::dot( vecLeft.t(), ( laplaceMatrix * vecRight ) );
  const double  result  =  0.5 * resultComplex.real(); // 0.5 to account for double-occupied orbitals due to spin
  return  result;
}


arma::sp_cx_dmat
SCFArmaSolver::getFMatrixKinetic_( const std::tuple<double, double>&  xRange,
                                   const std::tuple<double, double>&  yRange,
                                   const std::tuple<double, double>&  zRange
                                 ) const noexcept
{
  const arma::sp_cx_dmat  laplaceMatrix  =  getLaplaceMatrix_( xRange,
                                                               yRange,
                                                               zRange
                                                             );
  const unsigned&  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  const unsigned short * aPeriodicNumbers  =  geometry_->getPeriodicNumbers();

  std::list<std::size_t>  listIIndices;
  std::list<std::size_t>  listJIndices;
  std::list<double>  listMatrixElements;

  // to verify that the resulting sparse matrix will have the correct dimensions:
  bool  maxRowRegistered  =  false;
  bool  maxColRegistered  =  false;
  // < ^ >-product of each orbital of each atom
  // with each orbital of each other atom
  // and between orbitals within the same atom:
  std::size_t  iMatrixIndex  =  0;
  unsigned  iAtom1  =  0;
  for ( unsigned  iCoord1 = 0; iCoord1 < 3 * numAtoms; iCoord1 += 3 )
  {
    const double  xCenter1  =  aCoordinates[ iCoord1 ];
    const double  yCenter1  =  aCoordinates[ iCoord1 + 1 ];
    const double  zCenter1  =  aCoordinates[ iCoord1 + 2 ];
    const unsigned short  periodicNumber1  =  aPeriodicNumbers[ iAtom1 ];
    const std::size_t  numOrbitals1  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber1 ) );
    for ( std::size_t  iOrbital1 = 0; iOrbital1 < numOrbitals1; ++iOrbital1 )
    {
      std::size_t  jMatrixIndex  =  0;
      unsigned  iAtom2  =  0;
      for ( unsigned  iCoord2 = 0; iCoord2 < 3 * numAtoms; iCoord2 += 3 )
      {
        const double  xCenter2  =  aCoordinates[ iCoord2 ];
        const double  yCenter2  =  aCoordinates[ iCoord2 + 1 ];
        const double  zCenter2  =  aCoordinates[ iCoord2 + 2 ];
        const unsigned short  periodicNumber2  =  aPeriodicNumbers[ iAtom2 ];
        const std::size_t  numOrbitals2  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber2 ) );
        for ( std::size_t  iOrbital2 = 0; iOrbital2 < numOrbitals2; ++iOrbital2 )
        {
          const double  value  =  getFMatrixKineticValue_( laplaceMatrix,
                                                           xRange,
                                                           yRange,
                                                           zRange,
                                                           xCenter1,
                                                           yCenter1,
                                                           zCenter1,
                                                           xCenter2,
                                                           yCenter2,
                                                           zCenter2,
                                                           periodicNumber1,
                                                           periodicNumber2,
                                                           iOrbital1,
                                                           iOrbital2
                                                         );
          const double  eps  =  std::numeric_limits<double>::epsilon();
          if ( eps < value )
          {
            listIIndices.push_back( iMatrixIndex );
            listJIndices.push_back( jMatrixIndex );
            listMatrixElements.push_back( value );
            if ( numMolecularOrbitals_ - 1 == iMatrixIndex )
              maxRowRegistered  =  true;
            if ( numMolecularOrbitals_ - 1 == jMatrixIndex )
              maxColRegistered  =  true;
          }
          ++jMatrixIndex;
        } // for ( iOrbital2 )
        ++iAtom2;
      } // for ( iCoord2 )
      ++iMatrixIndex;
    } // for ( iOrbital1 )
    ++iAtom1;
  } // for ( iCoord1 )

  if ( false == maxRowRegistered )
  {
    listIIndices.push_back( numMolecularOrbitals_ - 1 );
    listJIndices.push_back( 0 );
    listMatrixElements.push_back( std::numeric_limits<double>::epsilon() );
  }
  if ( false == maxColRegistered )
  {
    listIIndices.push_back( 0 );
    listJIndices.push_back( numMolecularOrbitals_ - 1 );
    listMatrixElements.push_back( std::numeric_limits<double>::epsilon() );
  }
  // prepare data for the Armadillo sparce matrix constructor ( it there a better solution? ):
  const std::size_t  numNonZeroElementsMatrix  =  listIIndices.size();

  arma::cx_vec  vecMatrixValues( numNonZeroElementsMatrix, arma::fill::zeros );
  arma::umat  locationsInMatrix( 2, numNonZeroElementsMatrix, arma::fill::zeros );

  std::list<std::size_t>::iterator  itIIndices  =  listIIndices.begin();
  std::list<std::size_t>::iterator  itJIndices  =  listJIndices.begin();
  std::list<double>::iterator       itMatrixElements  =  listMatrixElements.begin();
  for ( std::size_t  i = 0;
        listIIndices.end() != itIIndices
     && listJIndices.end() != itJIndices
     && listMatrixElements.end() != itMatrixElements
     && i < numNonZeroElementsMatrix
        ;
        ++i
      )
  {
    locationsInMatrix( 0, i )  =  *itIIndices;
    locationsInMatrix( 1, i )  =  *itJIndices;
    const double  matrixElement  =  *itMatrixElements;
    vecMatrixValues( i )  =  arma::cx_double( matrixElement, 0 );
    ++itIIndices;
    ++itJIndices;
    ++itMatrixElements;
  }

  const arma::sp_cx_dmat  result( locationsInMatrix, vecMatrixValues );
  return  result;
}


arma::dvec
SCFArmaSolver::getFMatrixVecElectronsNucleiInteractionValues_( const std::tuple<double, double>&  xRange,
                                                               const std::tuple<double, double>&  yRange,
                                                               const std::tuple<double, double>&  zRange
                                                             ) const noexcept
{
  const unsigned&  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  const unsigned short * aPeriodicNumbers  =  geometry_->getPeriodicNumbers();

  arma::dvec  result( numMolecularOrbitals_, arma::fill::zeros );
  std::size_t  iMolecularOrbital  =  0;

  // account for interaction of each orbital of each atom
  // with each nucleus:
  for ( unsigned  iAtom1 = 0; iAtom1 < numAtoms; ++iAtom1 )
  {
    const unsigned short  periodicNumber1  =  aPeriodicNumbers[ iAtom1 ];
    const std::size_t  numOrbitals1  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber1 ) );
    for ( std::size_t  iOrbital1 = 0; iOrbital1 < numOrbitals1; ++iOrbital1 )
    {
      unsigned  iAtom2  =  0;
      double  interactionValue  =  0;
      for ( unsigned  iCoord2 = 0; iCoord2 < 3 * numAtoms; iCoord2 += 3 )
      {
        const double  xCenter2  =  aCoordinates[ iCoord2 ];
        const double  yCenter2  =  aCoordinates[ iCoord2 + 1 ];
        const double  zCenter2  =  aCoordinates[ iCoord2 + 2 ];
        const double  electronNucleusInteractionValue  =  getElectronNucleusInteractionValue_( xRange, yRange, zRange,
                                                                                               xCenter2, yCenter2, zCenter2,
                                                                                               periodicNumber1,
                                                                                               iOrbital1
                                                                                             );
        const unsigned short  periodicNumber2  =  aPeriodicNumbers[ iAtom2 ]; // == charge of the nucleus
        interactionValue  +=  ( periodicNumber2 * electronNucleusInteractionValue );
        ++iAtom2;
      }
      result( iMolecularOrbital )  =  interactionValue;
      ++iMolecularOrbital;
    }
  }
  return  result;
}


double
SCFArmaSolver::getElectronNucleusInteractionValue_( const std::tuple<double, double>&  xRange,
                                                    const std::tuple<double, double>&  yRange,
                                                    const std::tuple<double, double>&  zRange,
                                                    const double&                      xCenter1,
                                                    const double&                      yCenter1,
                                                    const double&                      zCenter1,
                                                    const unsigned short&              periodicNumber1,
                                                    const std::size_t&                 iOrbital1
                                                  ) const noexcept
{
  const double  margin  =  getMarginValue_();
  const double  lowX  =  std::get<0>( xRange ) - margin;
  const double  highX  =  std::get<1>( xRange ) + margin;
  const double  lowY  =  std::get<0>( yRange ) - margin;
  const double  highY  =  std::get<1>( yRange ) + margin;
  const double  lowZ  =  std::get<0>( zRange ) - margin;
  const double  highZ  =  std::get<1>( zRange ) + margin;

  Grid3D  grid( lowX, highX,
                lowY, highY,
                lowZ, highZ,
                resolution_,
                resolution_,
                resolution_
              );

  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;
  double  result  =  0;
  for ( std::size_t  iNode = 0; iNode < numNodes; ++iNode )
  {
    const std::size_t  iNodeZ  =  iNode / grid.numNodesX / grid.numNodesY;
    const std::size_t  iNodeY  =  ( iNode - iNodeZ * grid.numNodesX * grid.numNodesY ) / grid.numNodesX;
    const std::size_t  iNodeX  =  iNode - iNodeZ * grid.numNodesX * grid.numNodesY - iNodeY * grid.numNodesX;
    const double  x  =  lowX + iNodeX * resolution_;
    const double  y  =  lowY + iNodeY * resolution_;
    const double  z  =  lowZ + iNodeZ * resolution_;
    const double  val1  =  basisSet_->getValue( static_cast<short>( periodicNumber1 ),
                                                static_cast<short>( iOrbital1 ),
                                                xCenter1,
                                                yCenter1,
                                                zCenter1,
                                                x,
                                                y,
                                                z
                                              );
    const double  r2  =  ( x - xCenter1 ) * ( x - xCenter1 )
                       + ( y - yCenter1 ) * ( y - yCenter1 )
                       + ( z - zCenter1 ) * ( z - zCenter1 );
    const double  dist  =  std::sqrt( r2 );
    result  +=  ( val1 * val1 / dist );
  } // for ( iNode )
  return  result;
}


double
SCFArmaSolver::getElectronElectronInteractionValue_( const std::tuple<double, double>&  xRange,
                                                     const std::tuple<double, double>&  yRange,
                                                     const std::tuple<double, double>&  zRange,
                                                     const double&                      xCenter1,
                                                     const double&                      yCenter1,
                                                     const double&                      zCenter1,
                                                     const double&                      xCenter2,
                                                     const double&                      yCenter2,
                                                     const double&                      zCenter2,
                                                     const unsigned short&              periodicNumber1,
                                                     const std::size_t&                 iOrbital1,
                                                     const unsigned short&              periodicNumber2,
                                                     const std::size_t&                 iOrbital2
                                                   ) const noexcept
{
  const double  margin  =  getMarginValue_();
  const double  lowX  =  std::get<0>( xRange ) - margin;
  const double  highX  =  std::get<1>( xRange ) + margin;
  const double  lowY  =  std::get<0>( yRange ) - margin;
  const double  highY  =  std::get<1>( yRange ) + margin;
  const double  lowZ  =  std::get<0>( zRange ) - margin;
  const double  highZ  =  std::get<1>( zRange ) + margin;

  Grid3D  grid( lowX, highX,
                lowY, highY,
                lowZ, highZ,
                resolution_,
                resolution_,
                resolution_
              );

  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;
  double  result  =  0;
  for ( std::size_t  iNode = 0; iNode < numNodes; ++iNode )
  {
    const std::size_t  iNodeZ  =  iNode / grid.numNodesX / grid.numNodesY;
    const std::size_t  iNodeY  =  ( iNode - iNodeZ * grid.numNodesX * grid.numNodesY ) / grid.numNodesX;
    const std::size_t  iNodeX  =  iNode - iNodeZ * grid.numNodesX * grid.numNodesY - iNodeY * grid.numNodesX;
    const double  x1  =  lowX + iNodeX * resolution_;
    const double  y1  =  lowY + iNodeY * resolution_;
    const double  z1  =  lowZ + iNodeZ * resolution_;
    const double  val1  =  basisSet_->getValue( static_cast<short>( periodicNumber1 ),
                                                static_cast<short>( iOrbital1 ),
                                                xCenter1,
                                                yCenter1,
                                                zCenter1,
                                                x1,
                                                y1,
                                                z1
                                              );
    for ( std::size_t  iNode2 = 0; iNode2 < numNodes; ++iNode2 )
    {
      const std::size_t  iNodeZ  =  iNode2 / grid.numNodesX / grid.numNodesY;
      const std::size_t  iNodeY  =  ( iNode2 - iNodeZ * grid.numNodesX * grid.numNodesY ) / grid.numNodesX;
      const std::size_t  iNodeX  =  iNode2 - iNodeZ * grid.numNodesX * grid.numNodesY - iNodeY * grid.numNodesX;
      const double  x2  =  lowX + iNodeX * resolution_;
      const double  y2  =  lowY + iNodeY * resolution_;
      const double  z2  =  lowZ + iNodeZ * resolution_;

      const double  val2  =  basisSet_->getValue( static_cast<short>( periodicNumber2 ),
                                                  static_cast<short>( iOrbital2 ),
                                                  xCenter2,
                                                  yCenter2,
                                                  zCenter2,
                                                  x2,
                                                  y2,
                                                  z2
                                                );

      const double  r2  =  ( x1 - x2 ) * ( x1 - x2 )
                         + ( y1 - y2 ) * ( y1 - y2 )
                         + ( z1 - z2 ) * ( z1 - z2 );
      const double  dist  =  std::sqrt( r2 );
      result  +=  ( val1 * val1 * val2 * val2 / dist );
    } // for ( iNode2 )
  } // for ( iNode )
  return  result;
}


arma::sp_cx_dmat
SCFArmaSolver::getFMatrixElectronsElectronInteraction_( const std::tuple<double, double>&  xRange,
                                                        const std::tuple<double, double>&  yRange,
                                                        const std::tuple<double, double>&  zRange
                                                      ) const noexcept
{
  const unsigned&  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  const unsigned short * aPeriodicNumbers  =  geometry_->getPeriodicNumbers();

  arma::sp_cx_dmat  result;

  // account for interaction of each orbital of an atom
  // with each orbital of another atom:
  unsigned  iAtom1  =  0;
  std::size_t  iMolecularOrbital  =  0;
  for ( unsigned  iCoord1 = 0; iCoord1 < 3 * numAtoms; iCoord1 += 3 )
  {
    const double  xCenter1  =  aCoordinates[ iCoord1 ];
    const double  yCenter1  =  aCoordinates[ iCoord1 + 1 ];
    const double  zCenter1  =  aCoordinates[ iCoord1 + 2 ];
    const unsigned short  periodicNumber1  =  aPeriodicNumbers[ iAtom1 ];
    const std::size_t  numOrbitals1  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber1 ) );
    for ( std::size_t  iOrbital1 = 0; iOrbital1 < numOrbitals1; ++iOrbital1 )
    {
      unsigned  iAtom2  =  0;
      std::size_t  jMolecularOrbital  =  0;
      for ( unsigned  iCoord2 = 0; iCoord2 < 3 * numAtoms; iCoord2 += 3 )
      {
// <--- make 3 loops?!
        const double  xCenter2  =  aCoordinates[ iCoord2 ];
        const double  yCenter2  =  aCoordinates[ iCoord2 + 1 ];
        const double  zCenter2  =  aCoordinates[ iCoord2 + 2 ];
        const unsigned short  periodicNumber2  =  aPeriodicNumbers[ iAtom2 ];
        const std::size_t  numOrbitals2  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber2 ) );
        for ( std::size_t  iOrbital2 = 0; iOrbital2 < numOrbitals2; ++iOrbital2 )
        {
          if ( iMolecularOrbital == jMolecularOrbital )
            continue;
          const double  value  =  getElectronElectronInteractionValue_( xRange,
                                                                        yRange,
                                                                        zRange,
                                                                        xCenter1,
                                                                        yCenter1,
                                                                        zCenter1,
                                                                        xCenter2,
                                                                        yCenter2,
                                                                        zCenter2,
                                                                        periodicNumber1,
                                                                        iOrbital1,
                                                                        periodicNumber2,
                                                                        iOrbital2
                                                                      );
          const arma::cx_double  coefficient  =  mCoefficients_( iMolecularOrbital, jMolecularOrbital ); 
          result( iMolecularOrbital, jMolecularOrbital )  =  0.5 * std::conj( coefficient ) * coefficient * value; // 0.5 to avoid double-account of interactions
          ++jMolecularOrbital;
        }
        ++iAtom2;
      }
      ++iMolecularOrbital;
    } // for ( iOrbital1 )
    ++iAtom1;
  } // for ( iCoord1 )
  return  result;
}


arma::sp_cx_dmat
SCFArmaSolver::getFMatrix_( const std::tuple<double, double>&  xRange,
                            const std::tuple<double, double>&  yRange,
                            const std::tuple<double, double>&  zRange
                          ) const noexcept
{
  const arma::sp_cx_dmat  matrixKinetic  =  getFMatrixKinetic_( xRange, yRange, zRange );
  const arma::dvec  vecElectronsNucleiInteractions  =  getFMatrixVecElectronsNucleiInteractionValues_( xRange, yRange, zRange );
  const arma::sp_cx_dmat  matrixElectronElectron  =  getFMatrixElectronsElectronInteraction_( xRange, yRange, zRange );
  const arma::cx_dvec  vecExchangeCorrelation  =  getFMatrixExchangeCorrelationXAlpha_( xRange, yRange, zRange );
  arma::sp_cx_dmat  result  =  matrixKinetic + matrixElectronElectron;
  for ( unsigned i = 0; i < numMolecularOrbitals_; ++i )
  {
    result( i, i )  +=  vecElectronsNucleiInteractions( i );
    result( i, i )  +=  vecExchangeCorrelation( i );
  }
  return  result;
}


arma::cx_dvec
SCFArmaSolver::getFMatrixExchangeCorrelationXAlpha_( const std::tuple<double, double>&  xRange,
                                                     const std::tuple<double, double>&  yRange,
                                                     const std::tuple<double, double>&  zRange
                                                   ) const noexcept
{
  const unsigned&  numAtoms  =  geometry_->getNumAtoms();
  const double *  aCoordinates  =  geometry_->getCoordinates();
  const unsigned short * aPeriodicNumbers  =  geometry_->getPeriodicNumbers();

  arma::cx_dvec  result( numMolecularOrbitals_, arma::fill::zeros );

  for ( std::size_t  iMolecularOrbital = 0; iMolecularOrbital < numMolecularOrbitals_; ++iMolecularOrbital )
  {
    unsigned  iAtom1  =  0;
    std::size_t  jMolecularOrbital  =  0;
    for ( unsigned  iCoord1 = 0; iCoord1 < 3 * numAtoms; iCoord1 += 3 )
    {
      const double  xCenter1  =  aCoordinates[ iCoord1 ];
      const double  yCenter1  =  aCoordinates[ iCoord1 + 1 ];
      const double  zCenter1  =  aCoordinates[ iCoord1 + 2 ];
      const unsigned short  periodicNumber1  =  aPeriodicNumbers[ iAtom1 ];
      const std::size_t  numOrbitals1  =  basisSet_->getNumOrbitals( static_cast<short>( periodicNumber1 ) );
      for ( std::size_t  iOrbital1 = 0; iOrbital1 < numOrbitals1; ++iOrbital1 )
      {
        const double  value  =  getExchangeCorrelationXAlphaValue_( xRange,
                                                                    yRange,
                                                                    zRange,
                                                                    xCenter1,
                                                                    yCenter1,
                                                                    zCenter1,
                                                                    periodicNumber1,
                                                                    iOrbital1
                                                                  );
        const arma::cx_double  coefficient  =  mCoefficients_( iMolecularOrbital, jMolecularOrbital );
        constexpr double  alpha  =  ( 1. + 2. / 3 ) / 2; // see textbooks for empirical values	
        result( iMolecularOrbital )  +=  -9. / 8. * 3. / M_PI * alpha * std::conj( coefficient ) * coefficient * arma::cx_double( value, 0 );
        ++jMolecularOrbital;
      }
      ++iAtom1;
    }
  } // for ( iMolecularOrbital1 )
  return  result;
}


double
SCFArmaSolver::getExchangeCorrelationXAlphaValue_( const std::tuple<double, double>&  xRange,
                                                   const std::tuple<double, double>&  yRange,
                                                   const std::tuple<double, double>&  zRange,
                                                   const double&                      xCenter1,
                                                   const double&                      yCenter1,
                                                   const double&                      zCenter1,
                                                   const unsigned short&              periodicNumber1,
                                                   const std::size_t&                 iOrbital1
                                                 ) const noexcept
{
  const double  margin  =  getMarginValue_();
  const double  lowX  =  std::get<0>( xRange ) - margin;
  const double  highX  =  std::get<1>( xRange ) + margin;
  const double  lowY  =  std::get<0>( yRange ) - margin;
  const double  highY  =  std::get<1>( yRange ) + margin;
  const double  lowZ  =  std::get<0>( zRange ) - margin;
  const double  highZ  =  std::get<1>( zRange ) + margin;

  Grid3D  grid( lowX, highX,
                lowY, highY,
                lowZ, highZ,
                resolution_,
                resolution_,
                resolution_
              );

  const std::size_t  numNodes  =  grid.numNodesX * grid.numNodesY * grid.numNodesZ;
  double  result  =  0;
  for ( std::size_t  iNode = 0; iNode < numNodes; ++iNode )
  {
    const std::size_t  iNodeZ  =  iNode / grid.numNodesX / grid.numNodesY;
    const std::size_t  iNodeY  =  ( iNode - iNodeZ * grid.numNodesX * grid.numNodesY ) / grid.numNodesX;
    const std::size_t  iNodeX  =  iNode - iNodeZ * grid.numNodesX * grid.numNodesY - iNodeY * grid.numNodesX;
    const double  x1  =  lowX + iNodeX * resolution_;
    const double  y1  =  lowY + iNodeY * resolution_;
    const double  z1  =  lowZ + iNodeZ * resolution_;
    const double  val  =  basisSet_->getValue( static_cast<short>( periodicNumber1 ),
                                               static_cast<short>( iOrbital1 ),
                                               xCenter1,
                                               yCenter1,
                                               zCenter1,
                                               x1,
                                               y1,
                                               z1
                                             );
    result  +=  std::pow( val, 4. / 3. );
  } // for ( iNode )
  return  result;
}


} // namespace  cbc

