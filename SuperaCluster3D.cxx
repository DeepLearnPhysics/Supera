#ifndef __SUPERACLUSTER3D_CXX__
#define __SUPERACLUSTER3D_CXX__

#include "SuperaCluster3D.h"
#include "GenRandom.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaCluster3DProcessFactory __global_SuperaCluster3DProcessFactory__;

  SuperaCluster3D::SuperaCluster3D(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaCluster3D::configure(const PSet& cfg)
  {
		std::cerr << "Configuring SuperaCluster3D" << std::endl;
    SuperaBase::configure(cfg);
    _output_label   = cfg.get<std::string>("Cluster3DProducer");
  }

  void SuperaCluster3D::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaCluster3D::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    // Main cluster3d to be filled
    //auto& event_cluster_v = mgr.get_data<larcv::EventSparseTensor3D>(_output_label);
		auto& event2 = mgr.get_data<larcv::EventClusterVoxel3D>(_output_label);
		//auto const& meta = event_cluster_v.meta();
		auto const& meta = event2.meta();

    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();
		std::cout << meta.dump() << std::endl;

		larcv::VoxelSet cluster_spacepoint_v;


		auto const& spacepoint_v = LArData<supera::LArSpacePoint_t>();
		//cluster_spacepoint_v.resize(spacepoint_v.size()+1);
		LARCV_INFO() << "Processing SpacePoint array: " << spacepoint_v.size() << std::endl;
		for (size_t spacepoint_idx=0; spacepoint_idx < spacepoint_v.size(); ++spacepoint_idx) {
			//std::cout << spacepoint_idx << std::endl;
			auto const& spacepoint = spacepoint_v.at(spacepoint_idx);
			//std::cout << spacepoint << std::endl;
			auto coords = spacepoint.XYZ();
			auto error = spacepoint.ErrXYZ();
			//std::cout << error[0] << " " << error[1] << " " << error[2] << " " << error[3] << " " << error[4] << " " << error[5] << std::endl;
			//std::cout << coords << " " << spacepoint << std::endl;
			//std::cout << coords[0] << coords[1] << coords[2] << std::endl;
			//std::cout << meta.size_voxel_x() << meta.size_voxel_y() << meta.size_voxel_z() << std::endl;
			//std::cout << meta.min_x() << meta.max_x() << meta.min_y() << meta.max_y() << meta.min_z() << meta.max_z() << std::endl;
			//std::cout << (coords[0] - meta.min_x()) / meta.size_voxel_x() << std::endl;
			//std::cout << meta.id(coords[0], coords[1], coords[2]) << std::endl;
			larcv::VoxelID_t vox_id = meta.id(coords[0], coords[1], coords[2]);
			float vox_value = error[0];
			//std::cout << "SpacePoint " << vox_id << " " << coords[0] << std::endl;
			if (vox_id == larcv::kINVALID_VOXELID) {
				LARCV_DEBUG() << "Skipping SpacePoint from track id " << spacepoint.ID() << " pos=(" << coords[0] << "," << coords[1] << "," << coords[2] << ")" << std::endl;
				//std::cout << "hi" << std::endl;
				continue;
			}
			//cluster_spacepoint_v.push_back(larcv::Voxel(vox_id, 1.0));
			cluster_spacepoint_v.emplace(vox_id, vox_value, false);

		} // end for
		mgr.get_data<larcv::EventSparseTensor3D>("spacepoint").emplace(std::move(cluster_spacepoint_v), meta);

    return true;
  }

  void SuperaCluster3D::finalize()
  {}

}

#endif
