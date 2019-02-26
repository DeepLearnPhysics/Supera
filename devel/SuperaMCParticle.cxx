#ifndef __SUPERAMCPARTICLE_CXX__
#define __SUPERAMCPARTICLE_CXX__

#include "SuperaMCParticle.h"
#include "GenRandom.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaMCParticleProcessFactory __global_SuperaMCParticleProcessFactory__;

  SuperaMCParticle::SuperaMCParticle(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaMCParticle::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);

    _cluster3d_labels = cfg.get<std::vector<std::string> >("Cluster3DLabels");
    _tensor3d_labels  = cfg.get<std::vector<std::string> >("Tensor3DLabels");
 
  }

  void SuperaMCParticle::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaMCParticle::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto const& mcpart_v = LArData<supera::LArMCParticle_t>();
    

    // Register particle energy deposition coordinates
    auto const& sedep_v = LArData<supera::LArSimEnergyDeposit_t>();
    LARCV_INFO() << "Processing SimEnergyDeposit array: " << sedep_v.size() << std::endl;
    for(size_t sedep_idx=0; sedep_idx<sedep_v.size(); ++sedep_idx) {
      auto const& sedep = sedep_v.at(sedep_idx);
      larcv::Point3D pt;
      pt.x = sedep.X(); pt.y = sedep.Y(); pt.z=sedep.Z();
    }

    // Create Cluster3D
    for(auto const& name : _cluster3d_labels) {
      auto& cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(name);
      cluster3d.meta(meta);
    }
    // Create Tensor3D
    for(auto const& name : _tensor3d_labels) {
      auto& tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(name);
      tensor3d.meta(meta);
    }    
    return true;
  }
      
  void SuperaMCParticle::finalize()
  {}

}

#endif
