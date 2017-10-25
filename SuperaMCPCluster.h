/**
 * \file SuperaMCPCluster2D.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaMCPCluster2D
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMCPCLUSTER_H__
#define __SUPERAMCPCLUSTER_H__
#include "SuperaMCParticleTree.h"
#include "ParamsPixel2D.h"
#include "ImageMetaMaker.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaMCPCluster2D ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaMCPCluster2D : public supera::SuperaBase,
			     public supera::ImageMetaMaker,
			     public supera::ParamsPixel2D {
    
  public:
    
    /// Default constructor
    SuperaMCPCluster2D(const std::string name="SuperaMCPCluster2D");
    
    /// Default destructor
    ~SuperaMCPCluster2D(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::string _part_producer;
    MCParticleTree _mcpt;

  };

  /**
     \class larcv::SuperaMCPCluster2DFactory
     \brief A concrete factory class for larcv::SuperaMCPCluster2D
  */
  class SuperaMCPCluster2DProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCPCluster2DProcessFactory() { ProcessFactory::get().add_factory("SuperaMCPCluster2D",this); }
    /// dtor
    ~SuperaMCPCluster2DProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaMCPCluster2D(instance_name); }
  };

}
#endif
//#endif
//#endif
/** @} */ // end of doxygen group 

