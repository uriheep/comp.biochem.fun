
#include "grid3d.hpp"

#include <cmath>

namespace  cbc {

Grid3D::Grid3D() : numNodesX( 0 ),
                   numNodesY( 0 ),
                   numNodesZ( 0 ),
                   nodes( 0 )
{}


Grid3D::Grid3D( const double&  lowX,
                const double&  highX,
                const double&  lowY,
                const double&  highY,
                const double&  lowZ,
                const double&  highZ,
                const double&  stepX,
                const double&  stepY,
                const double&  stepZ
	      ) : numNodesX( 0 ),
                  numNodesY( 0 ),
                  numNodesZ( 0 ),
                  nodes( 0 )
{
  if ( lowX > highX
    || lowY > highY
    || lowZ > highZ
    || 0 >= stepX
    || 0 >= stepY
    || 0 >= stepZ
     )
  {
    return; // error
  }

  numNodesX  =  std::round( ( highX - lowX ) / stepX ) + 1;
  numNodesY  =  std::round( ( highY - lowY ) / stepY ) + 1;
  numNodesZ  =  std::round( ( highZ - lowZ ) / stepZ ) + 1;
}


} // namespace  cbc
