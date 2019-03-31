#ifndef __SUPERASEGSPACEPOINT_CXX__
#define __SUPERASEGSPACEPOINT_CXX__

#include "SuperaSegSpacePoint.h"
#include "GenRandom.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaSegSpacePointProcessFactory __global_SuperaSegSpacePointProcessFactory__;

  SuperaSegSpacePoint::SuperaSegSpacePoint(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaSegSpacePoint::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _output_label   = cfg.get<std::string>("Cluster3DProducer");
		_data_label = cfg.get<std::string>("DataProducer");
		_distance_threshold = cfg.get<float>("DistanceThreshold");
  }

  void SuperaSegSpacePoint::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaSegSpacePoint::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

		auto& event_spacepoint_v = mgr.get_data<larcv::EventSparseTensor3D>("spacepoint");
		auto& event_data_v = mgr.get_data<larcv::EventSparseTensor3D>(_data_label);
		auto const& meta = event_spacepoint_v.meta();

    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();

		larcv::VoxelSet cluster_spacepoint_v;
		//auto const& spacepoint_v = event_spacepoint_v.as_vector();
		//auto const& data_v = event_data_v.as_vector();
		for (auto const& spacepoint : event_spacepoint_v.as_vector()) {
			/*auto const& real_spacepoint = event_data_v.find(spacepoint.id());
			int vox_value = 1;
			if (real_spacepoint.id() == larcv::kINVALID_VOXELID) vox_value = 0;
			if (vox_value == 1) std::cout << spacepoint.id() << std::endl;*/
			int vox_value = 0;
			larcv::Point3D spacept = meta.position(spacepoint.id());
			std::vector<double> distances_v;
			double min_distance = 1e9;
			for (auto const& datapoint : event_data_v.as_vector()) {
				larcv::Point3D datapt = meta.position(datapoint.id());
				double d = spacept.distance(datapt);
				distances_v.push_back(d);
				if (d < min_distance) min_distance = d;
			}
			//std::cout << "Min distance " << min_distance << " " << event_data_v.as_vector().size() << std::endl;
			if (min_distance <= _distance_threshold) vox_value = 1;
			cluster_spacepoint_v.emplace(spacepoint.id(), vox_value, false);

		} // end for
		mgr.get_data<larcv::EventSparseTensor3D>("segspacepoint").emplace(std::move(cluster_spacepoint_v), meta);

    return true;
  }

  void SuperaSegSpacePoint::finalize()
  {}

}

#endif
