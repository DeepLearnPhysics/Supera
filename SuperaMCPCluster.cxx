#ifndef __SUPERAMCPCLUSTER_CXX__
#define __SUPERAMCPCLUSTER_CXX__

#include "SuperaMCPCluster.h"
#include "Voxel3DSlicer.h"
#include "PulledPork3DSlicer.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel2D.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"
#include "LAr2Image.h"

namespace larcv {

  static SuperaMCPClusterProcessFactory __global_SuperaMCPClusterProcessFactory__;

  SuperaMCPCluster::SuperaMCPCluster(const std::string name)
    : SuperaBase(name)
  {}

  void SuperaMCPCluster::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsPixel2D::configure(cfg);
    supera::ParamsVoxel3D::configure(cfg);
    //_mcpt.configure(cfg.get<supera::Config_t>("MCParticleTree"));
    _part_producer = cfg.get<std::string>("ParticleProducer");
    _target_projection = cfg.get<size_t>("TargetProjection",larcv::kINVALID_SIZE);
    _use_true_pos = cfg.get<bool>("UseTrue3DPosition");
  }

  void SuperaMCPCluster::initialize()
  { SuperaBase::initialize(); }

  bool SuperaMCPCluster::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }else if(supera::Voxel3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::Voxel3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());      
    }

    auto const& part_v = mgr.get_data<larcv::EventParticle>(_part_producer).as_vector();

    // create trackid=>clusterid mapping as a 1d array
    std::vector<size_t> trackid2cluster(1000, larcv::kINVALID_SIZE); // initialize to size 1000, cuz why not (cheaper than doing resize many times)
    for (auto const& part : part_v) {
      LARCV_INFO() << "Clustering particle (cluster id, track id, pdg) = (" 
		   << part.id() << "," << part.track_id() << "," << part.pdg_code() << ")" << std::endl;
      auto track_id = part.track_id();
      auto cluster_id = part.id();
      if (trackid2cluster.size() <= track_id) trackid2cluster.resize(track_id + 1, larcv::kINVALID_SIZE);
      trackid2cluster[track_id] = cluster_id;
    }

    // Register MCShower daughters if relevant
    auto const& mcshower_v = LArData<supera::LArMCShower_t>();
    for (auto const& mcs : mcshower_v) {
      if (mcs.TrackID() >= trackid2cluster.size()) continue;
      auto const cluster_id = trackid2cluster[mcs.TrackID()];
      if (cluster_id < 0) continue;
      for (auto const& daughter_id : mcs.DaughterTrackID()) {
        if (trackid2cluster.size() <= daughter_id) trackid2cluster.resize(daughter_id + 1, larcv::kINVALID_SIZE);
        trackid2cluster[daughter_id] = cluster_id;
      }
    }

    //
    // Is pixel2d cluster requested? shit
    //                                                                                                                                                                              
    if(!(OutPixel2DLabel().empty())) {
      auto meta_v = Meta();
      /*
      for (auto& meta : meta_v)
	meta.update(meta.rows() / RowCompressionFactor().at(meta.id()),
		    meta.cols() / ColCompressionFactor().at(meta.id()));
      */
      std::vector<larcv::ClusterPixel2D> clusters;
      for(auto const& meta : meta_v) {
	if(clusters.size() <= meta.id()) clusters.resize(meta.id()+1);

	clusters[meta.id()].resize(part_v.size()+1);
	clusters[meta.id()].meta(meta);
      }

      supera::SimCh2ClusterPixel2D(clusters, LArData<supera::LArSimCh_t>(), 
				   trackid2cluster, TimeOffset());

      auto& ev_pixel2d  = mgr.get_data<larcv::EventClusterPixel2D>(OutPixel2DLabel());

      for(size_t i=0; i<clusters.size(); ++i) {
	if( (int)(this->logger().level()) < (int)(larcv::msg::kINFO) ) {
	  auto const& plane_clusters = clusters[i];
	  LARCV_INFO() << "Plane meta: " << plane_clusters.meta().dump();
	  double total_sum=0;
	  for(size_t j=0; j<plane_clusters.as_vector().size(); ++j) {
	    auto const& cluster = plane_clusters.as_vector()[j];
	    double sum=0;
	    for(auto const& vox : cluster.as_vector())
	      sum += vox.value();
	    LARCV_INFO() << "Cluster " << j << " size " << cluster.as_vector().size() << " sum " << sum << std::endl;
	    total_sum += sum;
	  }
	  LARCV_INFO() << "Total sum: " << total_sum << std::endl;
	}
	ev_pixel2d.emplace(std::move(clusters[i]));
      }      
    }
    //
    // Is voxel3d cluster requested? dipshit
    //
    if(!(OutVoxel3DLabel().empty())) {
      auto meta3d = Meta3D();
      auto& ev_voxel3d = mgr.get_data<larcv::EventClusterVoxel3D>(OutVoxel3DLabel());
      ev_voxel3d.meta(meta3d);
      ev_voxel3d.resize(part_v.size()+1);
      supera::SimCh2ClusterVoxel3D(ev_voxel3d, LArData<supera::LArSimCh_t>(), 
				   trackid2cluster, TimeOffset(), _use_true_pos, _target_projection);
    }

    return true;
  }

  void SuperaMCPCluster::finalize()
  { SuperaBase::finalize(); }

}
#endif
