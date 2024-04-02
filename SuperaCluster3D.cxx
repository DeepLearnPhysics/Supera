#ifndef __SUPERACLUSTER3D_CXX__
#define __SUPERACLUSTER3D_CXX__

#include "SuperaCluster3D.h"
#include "GenRandom.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaCluster3DProcessFactory __global_SuperaCluster3DProcessFactory__;

  SuperaCluster3D::SuperaCluster3D(const std::string name) : SuperaBase(name) {}
  void SuperaCluster3D::configure(const PSet& cfg)
  {
    std::cerr << "Configuring SuperaCluster3D" << std::endl;
    SuperaBase::configure(cfg);
    _output_label = cfg.get<std::string>("Cluster3DProducer");
  }

  void SuperaCluster3D::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaCluster3D::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto& event2 = mgr.get_data<larcv::EventClusterVoxel3D>(_output_label);
    auto const& meta = event2.meta();

    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();

    larcv::VoxelSet cluster_spacepoint_v;

    auto const& spacepoint_v = LArData<supera::LArSpacePoint_t>();
    LARCV_INFO() << "Processing SpacePoint array: " << spacepoint_v.size() << std::endl;
    // Loop over space points and store them in cluster_spacepoint_v
    for (size_t spacepoint_idx = 0; spacepoint_idx < spacepoint_v.size(); ++spacepoint_idx) {
      auto const& spacepoint = spacepoint_v.at(spacepoint_idx);
      auto coords = spacepoint.XYZ();
      auto error = spacepoint.ErrXYZ();
      larcv::VoxelID_t vox_id = meta.id(coords[0], coords[1], coords[2]);
      float vox_value = error[0];
      if (vox_id == larcv::kINVALID_VOXELID) {
        LARCV_DEBUG() << "Skipping SpacePoint from id " << spacepoint.ID() << " pos=(" << coords[0]
                      << "," << coords[1] << "," << coords[2] << ")" << std::endl;
        continue;
      }
      cluster_spacepoint_v.emplace(vox_id, vox_value, false);

    } // end for
    // Write cluster_spacepoint_v to sparse3d_spacepoint_tree
    mgr.get_data<larcv::EventSparseTensor3D>("spacepoint")
      .emplace(std::move(cluster_spacepoint_v), meta);

    return true;
  }

  void SuperaCluster3D::finalize() {}

}

#endif
