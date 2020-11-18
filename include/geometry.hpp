#ifndef  CBC_GEOMETRY_HPP__
#define  CBC_GEOMETRY_HPP__

#define MAX_CHAR_ATOM_NAME 4
#define MAX_CHAR_RES_NAME 4

namespace  cbc {

class  Geometry {
  public:
    Geometry();
    ~Geometry();
    void  readCIF( const char * const  filename ) noexcept;
    // void  setCharges();
    // void  setMultiplicity();
    const unsigned&  getNumAtoms() const noexcept;
    const unsigned&  getNumChains() const noexcept;
    const unsigned&  getNumResidues() const noexcept;
    const unsigned short *  getPeriodicNumbers() const noexcept;
    double *  getCoordinates() const noexcept;
    double *  getCharges() const noexcept;
    const unsigned short *  getMultiplicity() const noexcept;
    unsigned short          getNumAtomsInResidue( const unsigned&  iResidue ) const noexcept;
    unsigned                getNumResiduesInChain( const unsigned&  iChain ) const noexcept;
    char                    getResidueName( const unsigned&  iResidue ) const noexcept;
    void  setCoordinates( const unsigned&  indexAtom,
                          const double&    x,
                          const double&    y,
                          const double&    z
                        ) noexcept;
  private:
    Geometry( const Geometry& );
    Geometry&  operator=( const Geometry& );
  private:
    short  getPeriodicNumber_( const char  (&atomName)[ MAX_CHAR_ATOM_NAME ] ) const noexcept;
    char   getResidueName_( const char  (&residueName)[ MAX_CHAR_RES_NAME ] ) const noexcept;
  private:
    // this design uses plain arrays intentionally:
    // for a smoother employment of CUDA in the future.
    unsigned          numAtoms_;
    unsigned          numChains_;
    unsigned          numResidues_;
    unsigned short  * aPeriodicNumbers_;
    double          * aCoordinates_;
    double          * aCharges_;
    unsigned short  * aMultiplicity_;
    unsigned short  * aNumAtomsInResidue_;
    unsigned        * aNumResiduesInChain_;
    char            * aResiduesNames_;
};

}

#endif
