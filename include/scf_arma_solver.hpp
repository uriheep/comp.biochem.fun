#ifndef  SCF_ARMADILLO_SOLVER_HPP__
#define  SCF_ARMADILLO_SOLVER_HPP__

#include <armadillo>

#include <list>
#include <tuple>


namespace  cbc {

class  Geometry;
class  BasisSetGTO;
class  Grid3D;

class  SCFArmaSolver {
  public:
    SCFArmaSolver( const Geometry * const     geometry,
                   const BasisSetGTO * const  basisSet
                 );
    ~SCFArmaSolver() noexcept;
    // resolution is in Angstroms:
    void  setResolution( const double&  delta ) noexcept;
    void  setMaxNumNodesInGrid( const std::size_t&  numNodes ) noexcept;
    const double&  getResolution() const noexcept;
    std::size_t    getNumNodesInGrid() const noexcept;
    void  run( const double&  tolerance ) noexcept;
    void  run( const std::size_t  numIterations ) noexcept;
    const double&  getEnergy() const noexcept;
  private:
    unsigned  getNumMolecularOrbitals_() const noexcept;
    std::list<unsigned>  getListNumAtomsOfEachType_() const noexcept;
    unsigned short  getNumBasisFunctions_() const noexcept;
    std::tuple<double, double>  getSpatialLimitsX_() const noexcept;
    std::tuple<double, double>  getSpatialLimitsY_() const noexcept;
    std::tuple<double, double>  getSpatialLimitsZ_() const noexcept;
    arma::sp_cx_dmat  getSMatrix_( const std::tuple<double, double>&  xRange,
                                   const std::tuple<double, double>&  yRange,
                                   const std::tuple<double, double>&  zRange
                                 ) const noexcept;
    double  getSMatrixValue_( const std::tuple<double, double>&  xRange,
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
                            ) const noexcept;
    arma::sp_cx_dmat  getLaplaceMatrix_( const std::tuple<double, double>&  xRange,
                                         const std::tuple<double, double>&  yRange,
                                         const std::tuple<double, double>&  zRange
                                       ) const noexcept;
    double  getFMatrixKineticValue_( const arma::sp_cx_dmat&            laplaceMatrix,
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
                                   ) const noexcept;
    arma::sp_cx_dmat  getFMatrixKinetic_( const std::tuple<double, double>&  xRange,
                                          const std::tuple<double, double>&  yRange,
                                          const std::tuple<double, double>&  zRange
                                        ) const noexcept;

    arma::dvec  getFMatrixVecElectronsNucleiInteractionValues_( const std::tuple<double, double>&  xRange,
                                                                const std::tuple<double, double>&  yRange,
                                                                const std::tuple<double, double>&  zRange
                                                              ) const noexcept;
    double  getElectronNucleusInteractionValue_( const std::tuple<double, double>&  xRange,
                                                 const std::tuple<double, double>&  yRange,
                                                 const std::tuple<double, double>&  zRange,
                                                 const double&                      xCenter1,
                                                 const double&                      yCenter1,
                                                 const double&                      zCenter1,
                                                 const unsigned short&              periodicNumber1,
                                                 const std::size_t&                 iOrbital1
                                               ) const noexcept;
    double  getElectronElectronInteractionValue_( const std::tuple<double, double>&  xRange,
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
                                                ) const noexcept;
    arma::cx_dmat  getFMatrixElectronsElectronInteraction_( const std::tuple<double, double>&  xRange,
                                                            const std::tuple<double, double>&  yRange,
                                                            const std::tuple<double, double>&  zRange
                                                          ) const noexcept;


    double  getFMatrixValue_( const std::tuple<double, double>&  xRange,
                              const std::tuple<double, double>&  yRange,
                              const std::tuple<double, double>&  zRange
                            ) const noexcept;
    void  initializeCoefficients_() noexcept;

  private:
    const Geometry        *geometry_;
    const BasisSetGTO     *basisSet_;
    unsigned              numMolecularOrbitals_;
    unsigned short        numBasisFunctions_;
    double                energy_;
    double                resolution_;
    std::size_t           maxNumNodesInGrid_;
    Grid3D                *aGrids_;
    //double                *aCoefficients_;
    arma::cx_dmat         mCoefficients_;
};

}

#endif
