
namespace  cbc {

struct  Grid3D {
  public:
    Grid3D();
    Grid3D( const double&  lowX,
            const double&  highX,
            const double&  lowY,
            const double&  highY,
            const double&  lowZ,
            const double&  highZ,
            const double&  stepX  =  1,
            const double&  stepY  =  1,
            const double&  stepZ  =  1
          );
  public:
    unsigned  numNodesX;
    unsigned  numNodesY;
    unsigned  numNodesZ;
    double    *nodes;
  private:
    // to avoid unintentional employment
    // of copy-ctors:
    Grid3D( const Grid3D& );
};

}
