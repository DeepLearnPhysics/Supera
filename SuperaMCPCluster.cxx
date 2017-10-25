#ifndef __SUPERAMCPCLUSTER_CXX__
#define __SUPERAMCPCLUSTER_CXX__

#include "SuperaMCPCluster.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventPixel2D.h"
#include "LAr2Image.h"

namespace larcv {

  static SuperaMCPClusterProcessFactory __global_SuperaMCPClusterProcessFactory__;

  SuperaMCPCluster::SuperaMCPCluster(const std::string name)
    : SuperaMCROI(name)
  {}

  void SuperaMCPCluster::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsPixel2D::configure(cfg);
    _mcpt.configure(cfg.get<supera::Config_t>("MCParticleTree"));
    _part_producer = cfg.get<std::string>("ParticleProducer");
  }

  void SuperaMCPCluster::initialize()
  { SuperaMCROI::initialize(); }

  bool SuperaMCPCluster::process(IOManager& mgr)
  {
    SuperaMCROI::process(mgr);

    auto& ev_pixel2d   = mgr.get_data<larcv::EventPixel2D>(OutPixel2DLabel());
    auto const& part_v = mgr.get_data<larcv::EventParticle>(_part_producer).particle_array();
    
    // create trackid=>clusterid mapping as a 1d array
    std::vector<int> trackid2cluster(1000,-1); // initialize to size 1000, cuz why not (cheaper than doing resize many times)
    for(auto const& part : part_v) {
      auto track_id = part.track_id();
      auto cluster_id = part.id();
      if(trackid2cluster.size() <= track_id) trackid2cluster.resize(track_id+1,-1);
      trackid2cluster[track_id] = cluster_id;
    }

    // Register MCShower daughters if relevant
    auto const& mcshower_v = LArData<supera::LArMCShower_t>();
    for(auto const& mcs : mcshower_v) {
      if(mcs.TrackID() >= trackid2cluster.size()) continue;
      auto const cluster_id = trackid2cluster[mcs.TrackID()];
      if(cluster_id<0) continue;
      for(auto const& daughter_id : mcs.Daughters()) {
	if(trackid2cluster.size() <= daughter_id) trackid2cluster.resize(daughter_id+1,-1);
	trackid2cluster[daughter_id] = cluster_id;
      }
    }

    auto meta_v = Meta();
    for(auto& meta : meta_v) 
      meta.update(meta.rows() / RowCompressionFactor().at(meta.plane()),
		  meta.cols() / ColCompressionFactor().at(meta.plane()));

    auto res = supera::SimCh2Pixel2DCluster(meta_v,LArData<supera::LArSimCh_t>(),TimeOffset());

    ev_pixel2d.emplace(std::move(res));

    return true;
  }

  void SuperaMCPCluster::finalize()
  { SuperaMCROI::finalize(); }

}
#endif
