#ifndef  __BASIS_SET_STO_nG_HPP__
#define  __BASIS_SET_STO_nG_HPP__

#include <cstddef>
#include <string>


namespace  cbc {

class  BasisSetGTO {
  public:
    BasisSetGTO( const std::string&  filename );
    ~BasisSetGTO();
    const std::size_t&  getNumChemElements() const noexcept;
    std::size_t         getNumOrbitals( const short& ) const noexcept;
    double              getValue( const short&   indexElement,
                                  const short&   indexOrbital,
                                  const double&  xCenter,
                                  const double&  yCenter,
                                  const double&  zCenter,
                                  const double&  x,
                                  const double&  y,
                                  const double&  z
                                ) const noexcept;
  private:
    BasisSetGTO( const BasisSetGTO& );
    double  getSOrbitalPrimitive_( const short&   iElement,
                                   const short&   iOrbital,
                                   const short&   iPrimitive,
                                   const double&  r2
                                 ) const noexcept;
    double  getPOrbitalPrimitive_( const short&   iElement,
                                   const short&   iOrbital,
                                   const short&   iPrimitive,
                                   const double&  r2,
                                   const double&  x
                                 ) const noexcept;
    double  getDOrbitalPrimitive_( const short&   iElement,
                                   const short&   iOrbital,
                                   const short&   iPrimitive,
                                   const double&  r2
                                 ) const noexcept;
  private:
    struct  Orbital_ {
      public:
        Orbital_() : numPrimitives( -1 ),
                     angularMomentum( -1 ),
                     arrExponents( nullptr ),
                     arrCoeffs( nullptr )
        { }
      public:
        short   numPrimitives;
        short   angularMomentum;
        double  *arrExponents;
        double  *arrCoeffs;
    };
    struct  Element_ {
      public:
        Element_() : iPeriodic( -1 ),
                     numOrbitals( 1000 ),
                     arrOrbitals( nullptr )
        { }
      public:
        short        iPeriodic;
        std::size_t  numOrbitals;
        Orbital_     *arrOrbitals;
    };
  private:
    std::size_t  numChemElements_;
    Element_ *   arrElements_;
};

}

#endif
