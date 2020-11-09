// rapidjson/example/simpledom/simpledom.cpp`

#include "basis_set_sto-nG.hpp"

#include "rapidjson/document.h"
#include <cstdio>
#include <cmath>


namespace  cbc {

BasisSetGTO::BasisSetGTO( const std::string&  filename ) : numChemElements_( 0 ),
                                                           arrElements_( nullptr )
{
  if ( true == filename.empty() )
    return; // error

  std::FILE  *pFile  =  nullptr;
  pFile  =  std::fopen( filename.c_str(), "r" );
  if ( nullptr == pFile )
    printf( "BasisSetGTO: ctor error #0: %s\n", filename.c_str() );
  if ( 0 != std::fseek( pFile, 0, SEEK_END ) )
    printf( "BasisSetGTO: ctor error #1\n" );
  const std::size_t  fsize  =  std::ftell( pFile );
  if ( 0 != std::fseek( pFile, 0, SEEK_SET ) )
    printf( "BasisSetGTO: ctor error #2\n" );

  char * string  =  new char [ fsize + 1 ];
  std::fread( string, fsize, 1, pFile );
  std::fclose( pFile );
  pFile  =  nullptr;

  // parse a JSON string into DOM.
  rapidjson::Document  doc;
  doc.Parse( string );

  delete [] string;
  string  =  nullptr;

  // ******
  if ( false == doc.IsObject() )
    return; // error

  const rapidjson::Value&  elems  =  doc[ "elements" ];
    
  if ( false == elems.IsObject() )
    return; // error

  for ( rapidjson::Value::ConstMemberIterator  itElem = elems.MemberBegin();
        elems.MemberEnd() != itElem;
        ++itElem
      )
  {
    ++numChemElements_;
  }

  arrElements_  =  new Element_ [ numChemElements_ ];
  short  indexElem  =  0;

  for ( rapidjson::Value::ConstMemberIterator  itElem = elems.MemberBegin();
        elems.MemberEnd() != itElem;
        ++itElem
      )
  {
    const rapidjson::Value&  key  =  itElem->name;
    arrElements_[ indexElem ].iPeriodic  =  std::atoi( key.GetString() );

    arrElements_[ indexElem ].numOrbitals  =  0;

    const rapidjson::Value&  val  =  itElem->value;
    const rapidjson::Value&  elShells  =  val[ "electron_shells" ];

    if ( true == elShells.IsArray() )
    {
      for ( rapidjson::Value::ConstValueIterator  itShell = elShells.Begin();
            elShells.End() != itShell;
            ++itShell
          )
      {
        const rapidjson::Value&  arrAngularMomentum  =  (*itShell)[ "angular_momentum" ];
        arrElements_[ indexElem ].numOrbitals  +=  arrAngularMomentum.Size();
      }

      const std::size_t  numOrbitals  =  arrElements_[ indexElem ].numOrbitals;

      arrElements_[ indexElem ].arrOrbitals  =  new Orbital_ [ numOrbitals ];

      std::size_t  iOrbital  =  0;
      for ( rapidjson::Value::ConstValueIterator  itShell = elShells.Begin();
            elShells.End() != itShell;
            ++itShell
          )
      {
        const rapidjson::Value&  arrExp  =  (*itShell)[ "exponents" ];
        const rapidjson::Value&  arrCoeffs  =  (*itShell)[ "coefficients" ];
        if ( true == arrExp.IsArray()
          && true == arrCoeffs.IsArray()
           )
        {
          const short  numPrimitives  =  arrExp.Size();
          arrElements_[ indexElem ].arrOrbitals[ iOrbital ].numPrimitives  =  numPrimitives;
          // check it here due to the format of input JSON files
          // where exponents are given only once for all the
          // s, p and d orbitals:
          if ( nullptr == arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrExponents )
          {
            arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrExponents  =  new double [ numPrimitives ];
            arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrCoeffs  =  new double [ numPrimitives ];
          }
          // exponents:
          for ( rapidjson::SizeType i = 0; i < arrExp.Size(); ++i )
          {
            if ( true == arrExp[ i ].IsString()
              && 0 < std::strlen( arrExp[ i ].GetString() )
               )
            {
              const double  exp  =  std::atof( arrExp[ i ].GetString() );
              arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrExponents[ i ]  =  exp;
            }
          }
          const rapidjson::Value&  arrAngularMomentum  =  (*itShell)[ "angular_momentum" ];
          // coefficients:
          for ( rapidjson::SizeType i = 0; i < arrCoeffs.Size(); ++i )
          {
//if ( 1 == ang...
            if ( true == arrAngularMomentum[ i ].IsString()
              && 0 < std::strlen( arrAngularMomentum[ i ].GetString() )
               )
            {
              const short  angMomentum  =  std::atoi( arrAngularMomentum[ i ].GetString() );
              arrElements_[ indexElem ].arrOrbitals[ iOrbital ].angularMomentum  =  angMomentum;
            }
            const rapidjson::Value&  arr  =  arrCoeffs[ i ];
            for ( rapidjson::SizeType j = 0; j < arr.Size(); ++j )
            {
              if ( true == arr[ j ].IsString()
                && 0 < std::strlen( arr[ j ].GetString() )
                 )
              {
                const double  coeff  =  std::atof( arr[ j ].GetString() );
                arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrCoeffs[ j ]  =  coeff;
              }
            }
            ++iOrbital;
            if ( 1 < arrCoeffs.Size()
              && iOrbital < arrElements_[ indexElem ].numOrbitals
               )
            {
              if ( nullptr == arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrExponents )
              {
                arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrExponents  =  new double [ numPrimitives ];
                arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrCoeffs  =  new double [ numPrimitives ];
              }
              // fill the missing exponents on a new orbital:
              for ( short j = 0; j < numPrimitives; ++j )
                arrElements_[ indexElem ].arrOrbitals[ iOrbital ].arrExponents[ j ]  =  arrElements_[ indexElem ].arrOrbitals[ iOrbital - 1 ].arrExponents[ j ];
            }
          }
        } // if ( true == IsArray() )
      } // for ( itShell )
      ++indexElem;
    } // if ( IsArray() )
  } // for ( itElem )
}

BasisSetGTO::~BasisSetGTO()
{
  for ( std::size_t iElem = 0; iElem < numChemElements_; ++iElem )
  {
    const std::size_t  numOrbitals  =  arrElements_[ iElem ].numOrbitals;
    for ( std::size_t iOrbit = 0; iOrbit < numOrbitals; ++iOrbit )
    {
      delete [] arrElements_[ iElem ].arrOrbitals[ iOrbit ].arrExponents;
      arrElements_[ iElem ].arrOrbitals[ iOrbit ].arrExponents  =  nullptr;
      delete [] arrElements_[ iElem ].arrOrbitals[ iOrbit ].arrCoeffs;
      arrElements_[ iElem ].arrOrbitals[ iOrbit ].arrCoeffs  =  nullptr;
    }
    delete [] arrElements_[ iElem ].arrOrbitals;
    arrElements_[ iElem ].arrOrbitals  =  nullptr;
  }
  delete [] arrElements_;
  arrElements_  =  nullptr;
}

const std::size_t&
BasisSetGTO::getNumChemElements() const noexcept { return  numChemElements_; }

std::size_t
BasisSetGTO::getNumOrbitals( const short&  iElement ) const noexcept
{ 
  if ( nullptr == arrElements_
    || 0 > iElement
    || static_cast<std::size_t>( iElement ) >= numChemElements_
     )
  {
    return  0; // error
  }

  return  arrElements_[ iElement ].numOrbitals;
}

double
BasisSetGTO::getValue( const short&   indexElement,
                       const short&   indexOrbital,
                       const double&  xCenter,
                       const double&  yCenter,
                       const double&  zCenter,
                       const double&  x,
                       const double&  y,
                       const double&  z
                     ) const noexcept
{
  if ( 0 > indexElement
    || 0 > indexOrbital
    || static_cast<std::size_t>( indexElement ) >= numChemElements_
    || static_cast<std::size_t>( indexOrbital ) >= arrElements_[ indexElement ].numOrbitals
     )
  {
    return  0; // error
  }

  const short&  numPrimitives  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].numPrimitives;
  const double  r2  =  ( x - xCenter ) * ( x - xCenter )
                     + ( y - yCenter ) * ( y - yCenter )
                     + ( z - zCenter ) * ( z - zCenter );
//  const double  r  =  std::sqrt( r2 );
  const short&   angMomentum  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].angularMomentum;
  double  result  =  0;
  for ( short iPrim = 0; iPrim < numPrimitives; ++iPrim )
  {
    if ( 0 == angMomentum )
      result  +=  getSOrbitalPrimitive_( indexElement, indexOrbital, iPrim, r2 );
    if ( 1 == angMomentum )
    {
      result  +=  getPOrbitalPrimitive_( indexElement, indexOrbital, iPrim, r2, x - xCenter );
//which of 3 to return?
      result  +=  getPOrbitalPrimitive_( indexElement, indexOrbital, iPrim, r2, y - yCenter );
      result  +=  getPOrbitalPrimitive_( indexElement, indexOrbital, iPrim, r2, z - zCenter );
    }
    if ( 2 == angMomentum )
    {
      result  +=  getDOrbitalPrimitive_( indexElement, indexOrbital, iPrim, r2 );
    }
  }
  return  result;
}


double
BasisSetGTO::getSOrbitalPrimitive_( const short&   indexElement,
                                    const short&   indexOrbital,
                                    const short&   iPrim,
                                    const double&  r2
                                  ) const noexcept
{
  const double  coeff  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].arrCoeffs[ iPrim ];
  const double  exponent  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].arrExponents[ iPrim ];
  return  coeff * std::exp( -exponent * r2 );
}

double
BasisSetGTO::getPOrbitalPrimitive_( const short&   indexElement,
                                    const short&   indexOrbital,
                                    const short&   iPrim,
                                    const double&  r2,
                                    const double&  delta
                                  ) const noexcept
{
  const double  coeff  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].arrCoeffs[ iPrim ];
  const double  exponent  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].arrExponents[ iPrim ];
  return  coeff * delta * std::exp( -exponent * r2 );
}

double
BasisSetGTO::getDOrbitalPrimitive_( const short&   indexElement,
                                    const short&   indexOrbital,
                                    const short&   iPrim,
                                    const double&  r2
                                  ) const noexcept
{
  const double  coeff  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].arrCoeffs[ iPrim ];
  const double  exponent  =  arrElements_[ indexElement ].arrOrbitals[ indexOrbital ].arrExponents[ iPrim ];
  return  coeff * r2 * std::exp( -exponent * r2 );
}

} // namespace  cbc
