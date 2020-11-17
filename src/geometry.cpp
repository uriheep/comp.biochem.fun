

#include "geometry.hpp"

#include "readcif.h"
#include <iostream>
#include <limits>
#include <cstring>

// declarations of auxiliary functions
// borrowed from the READCIF library:
// file  'readcif_example.cpp'
// ( to be used in method  Geometry::readCIF() )

using std::string;
using std::vector;
using namespace readcif;

#define MAX_CHAR_CHAIN_ID 4

static const bool Required = true;

struct Atom
{
    Atom() {
        clear();
    }
    void clear() {
        memset(atom_name, 0, MAX_CHAR_ATOM_NAME);
        memset(residue_name, 0, MAX_CHAR_RES_NAME);
        memset(chain_id, 0, MAX_CHAR_CHAIN_ID);
    }
    char element;
    char atom_name[MAX_CHAR_ATOM_NAME];
    char residue_name[MAX_CHAR_RES_NAME];
    char chain_id[MAX_CHAR_CHAIN_ID];
    int residue_num;
    float x, y, z;
};

struct ExtractCIF: CIFFile {
    ExtractCIF();
    void parse_atom_site();
    std::vector<Atom> atoms;
};

/* **************************** */


namespace  cbc {


Geometry::Geometry() : numAtoms_( 0 ),
                       numChains_( 0 ),
                       numResidues_( 0 ),
                       aPeriodicNumbers_( nullptr ),
                       aCoordinates_( nullptr ),
                       aCharges_( nullptr ),
                       aMultiplicity_( nullptr ),
                       aNumAtomsInResidue_( nullptr ),
                       aNumResiduesInChain_( nullptr ),
                       aResiduesNames_( nullptr )
{ }

void
Geometry::readCIF( const char * const  filename ) noexcept
{
  if ( nullptr == filename )
    return;

  ExtractCIF  extract;

  try {
    extract.parse_file( filename );
  }
  catch ( std::exception&  e) {
    std::cerr << e.what() << '\n';
    return;
  }

  numAtoms_  =  extract.atoms.size();
  numResidues_  =  0;
  numChains_  =  0;

  delete [] aPeriodicNumbers_;
  delete [] aCoordinates_;
  delete [] aCharges_;
  delete [] aMultiplicity_;
  delete [] aNumAtomsInResidue_;
  delete [] aNumResiduesInChain_;
  delete [] aResiduesNames_;

  aPeriodicNumbers_  =  new unsigned short [ numAtoms_ ];
  aCoordinates_  =  new double [ 3 * numAtoms_ ];

  int  residueNumTmp  =  std::numeric_limits<int>::min(); // initial
  char  chainIDTmp[ MAX_CHAR_CHAIN_ID ]  =  { '0', '0', '0', '\0' }; // initial
  for ( std::size_t  iAtom = 0; iAtom < numAtoms_; ++iAtom )
  {
    const Atom  &atom  =  extract.atoms[ iAtom ];
    aPeriodicNumbers_[ iAtom ]  =  getPeriodicNumber_( atom.atom_name );
    aCoordinates_[ 3 * iAtom ]      =  atom.x;
    aCoordinates_[ 3 * iAtom + 1 ]  =  atom.y;
    aCoordinates_[ 3 * iAtom + 2 ]  =  atom.z;
    const int  residueNum  =  atom.residue_num;
    if ( residueNum != residueNumTmp )
    {
      ++numResidues_;
      residueNumTmp  =  residueNum;
    }
    if ( 0 != std::strcmp( chainIDTmp, atom.chain_id ) )
    {
      ++numChains_;
      std::strcpy( chainIDTmp, atom.chain_id );
    }
  }
  aNumAtomsInResidue_   =  new unsigned short [ numResidues_ ];
  aResiduesNames_       =  new char [ numResidues_ ];
  aNumResiduesInChain_  =  new unsigned short [ numChains_ ];

  residueNumTmp  =  extract.atoms[ 0 ].residue_num; // initial
  std::strcmp( chainIDTmp, extract.atoms[ 0 ].chain_id ); // initial
  unsigned short  numAtomsInResidue   =  0;
  unsigned short  numResiduesInChain  =  0;
  std::size_t     iResidue  =  0;
  std::size_t     iChain    =  0;
  for ( std::size_t  iAtom = 0; iAtom < numAtoms_; ++iAtom )
  {
    const Atom  &atom  =  extract.atoms[ iAtom ];
    ++numAtomsInResidue;
    const int  residueNum  =  atom.residue_num;
    if ( residueNum != residueNumTmp )
    {
      aNumAtomsInResidue_[ iResidue ]  =  numAtomsInResidue;
      aResiduesNames_[ iResidue ]  =  getResidueName_( atom.residue_name );
      ++iResidue;
      numAtomsInResidue  =  0;
      residueNumTmp  =  residueNum;
    }
    ++numResiduesInChain;
    if ( 0 != std::strcmp( chainIDTmp, atom.chain_id ) )
    {
      aNumResiduesInChain_[ iChain ]  =  numResiduesInChain;
      ++iChain;
      numResiduesInChain  =  0;
      std::strcpy( chainIDTmp, atom.chain_id );
    }
  }
}

Geometry::~Geometry() noexcept
{
  delete [] aPeriodicNumbers_;
  aPeriodicNumbers_  =  nullptr;
  delete [] aCoordinates_;
  aCoordinates_  =  nullptr;
  delete [] aCharges_;
  aCharges_  =  nullptr;
  delete [] aMultiplicity_;
  aMultiplicity_  =  nullptr;
  delete [] aNumAtomsInResidue_;
  aNumAtomsInResidue_  =  nullptr;
  delete [] aNumResiduesInChain_;
  aNumResiduesInChain_  =  nullptr;
  delete [] aResiduesNames_;
  aResiduesNames_  =  nullptr;
}

unsigned
Geometry::getNumAtoms() const noexcept { return  numAtoms_; }

const unsigned short *
Geometry::getPeriodicNumbers() const noexcept { return  aPeriodicNumbers_; }

double *
Geometry::getCoordinates() const noexcept { return  aCoordinates_; }

double *
Geometry::getCharges() const noexcept { return  aCharges_; }

const unsigned short *
Geometry::getMultiplicity() const noexcept { return  aMultiplicity_; }

void
Geometry::setCoordinates( const unsigned&  indexAtom,
                          const double&    x,
                          const double&    y,
                          const double&    z
                        ) noexcept
{
  if ( nullptr == aCoordinates_
    || indexAtom >= numAtoms_
     )
    return;

  aCoordinates_[ 3 * indexAtom ]      =  x;
  aCoordinates_[ 3 * indexAtom + 1 ]  =  y;
  aCoordinates_[ 3 * indexAtom + 2 ]  =  z;
}

short
Geometry::getPeriodicNumber_( const char  (&atomName)[ MAX_CHAR_ATOM_NAME ] ) const noexcept
{
  if ( 'H' == atomName[ 0 ] )
    return  1;
  if ( 'C' == atomName[ 0 ] )
    return  6;
  if ( 'N' == atomName[ 0 ] )
    return  7;
  if ( 'O' == atomName[ 0 ] )
    return  8;
  if ( 'S' == atomName[ 0 ] )
    return  16;
  return  0; // error code
}

char
Geometry::getResidueName_( const char  (&residueName)[ MAX_CHAR_RES_NAME ] ) const noexcept
{
  if ( 0 == std::strcmp( "ALA\0", residueName ) )  return  'A';
  if ( 0 == std::strcmp( "ARG\0", residueName ) )  return  'R';
  if ( 0 == std::strcmp( "ASN\0", residueName ) )  return  'N';
  if ( 0 == std::strcmp( "ASP\0", residueName ) )  return  'D';
  if ( 0 == std::strcmp( "CYS\0", residueName ) )  return  'C';
  if ( 0 == std::strcmp( "GLN\0", residueName ) )  return  'Q';
  if ( 0 == std::strcmp( "GLU\0", residueName ) )  return  'E';
  if ( 0 == std::strcmp( "GLY\0", residueName ) )  return  'G';
  if ( 0 == std::strcmp( "HIS\0", residueName ) )  return  'H';
  if ( 0 == std::strcmp( "ILE\0", residueName ) )  return  'I';
  if ( 0 == std::strcmp( "LEU\0", residueName ) )  return  'L';
  if ( 0 == std::strcmp( "LYS\0", residueName ) )  return  'K';
  if ( 0 == std::strcmp( "MET\0", residueName ) )  return  'M';
  if ( 0 == std::strcmp( "PHE\0", residueName ) )  return  'F';
  if ( 0 == std::strcmp( "PRO\0", residueName ) )  return  'P';
  if ( 0 == std::strcmp( "SER\0", residueName ) )  return  'S';
  if ( 0 == std::strcmp( "THR\0", residueName ) )  return  'T';
  if ( 0 == std::strcmp( "TRP\0", residueName ) )  return  'W';
  if ( 0 == std::strcmp( "TYR\0", residueName ) )  return  'Y';
  if ( 0 == std::strcmp( "VAL\0", residueName ) )  return  'V';
  if ( 0 == std::strcmp( "ASX\0", residueName ) )  return  'B';
  if ( 0 == std::strcmp( "GLX\0", residueName ) )  return  'Z';
  return  '0'; // cannot identify residue type
}

} // namespace cbc


ExtractCIF::ExtractCIF()
{
    register_heuristic_stylized_detection();
#if 0
    using std::placeholder;
    register_category("atom_site", 
                      std::bind(&ExtractCIF::parse_atom_site, this, _1));
#else
    // Personal preference, I like lambda functions better.
    // The lambda functions are needed because parse_XXXX
    // are member functions.
    register_category("atom_site", 
                      [this] () {
                          parse_atom_site();
                      });
#endif
}

void
ExtractCIF::parse_atom_site()
{
    CIFFile::ParseValues pv;
    pv.reserve(10);
    Atom atom;
    pv.emplace_back(get_column("type_symbol", Required),
                    [&atom] (const char* start) {
                        atom.element = *start;
                    });
    pv.emplace_back(get_column("label_atom_id", Required),
                    [&atom] (const char* start, const char* end) {
                        size_t count = end - start;
                        if (count > MAX_CHAR_ATOM_NAME)
                            count = MAX_CHAR_ATOM_NAME;
                        strncpy(atom.atom_name, start, count);
                    });
    pv.emplace_back(get_column("label_comp_id", Required),
                    [&atom] (const char* start, const char* end) {
                        size_t count = end - start;
                        if (count > MAX_CHAR_RES_NAME)
                            count = MAX_CHAR_RES_NAME;
                        strncpy(atom.residue_name, start, count);
                    });
    pv.emplace_back(get_column("label_asym_id"),
                    [&atom] (const char* start, const char* end) {
                        size_t count = end - start;
                        if (count > MAX_CHAR_CHAIN_ID)
                            count = MAX_CHAR_CHAIN_ID;
                        strncpy(atom.chain_id, start, count);
                    });
    pv.emplace_back(get_column("label_seq_id", Required),
                    [&atom] (const char* start) {
                        atom.residue_num = readcif::str_to_int(start);
                    });
    // x, y, z are not required by mmCIF, but are by us
    pv.emplace_back(get_column("Cartn_x", Required),
                    [&atom] (const char* start) {
                        atom.x = readcif::str_to_float(start);
                    });
    pv.emplace_back(get_column("Cartn_y", Required),
                    [&atom] (const char* start) {
                        atom.y = readcif::str_to_float(start);
                    });
    pv.emplace_back(get_column("Cartn_z", Required),
                    [&atom] (const char* start) {
                        atom.z = readcif::str_to_float(start);
                    });
    while (parse_row(pv)) {
        atoms.push_back(atom);
        atom.clear();
    }
}

