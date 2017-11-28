/**
 * \file SuperaMCPCluster.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaMCPCluster
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMCPCLUSTER_H__
#define __SUPERAMCPCLUSTER_H__
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "MCParticleTree.h"
#include "ParamsPixel2D.h"
#include "ParamsVoxel3D.h"
#include "ImageMetaMaker.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaMCPCluster ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaMCPCluster : public SuperaBase,
			   public supera::ImageMetaMaker,
			   public supera::ParamsPixel2D,
			   public supera::ParamsVoxel3D {
    
  public:
    
    /// Default constructor
    SuperaMCPCluster(const std::string name="SuperaMCPCluster");
    
    /// Default destructor
    ~SuperaMCPCluster(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::string _part_producer;
    ProjectionID_t _target_projection;
    supera::MCParticleTree _mcpt;
    bool _use_true_pos;
  };

  /**
     \class larcv::SuperaMCPClusterFactory
     \brief A concrete factory class for larcv::SuperaMCPCluster
  */
  class SuperaMCPClusterProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCPClusterProcessFactory() { ProcessFactory::get().add_factory("SuperaMCPCluster",this); }
    /// dtor
    ~SuperaMCPClusterProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaMCPCluster(instance_name); }
  };

}
#endif
//#endif
//#endif
/** @} */ // end of doxygen group 

