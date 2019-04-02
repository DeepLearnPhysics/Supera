#ifndef __SUPERASIMENERGYDEPOSIT_CXX__
#define __SUPERASIMENERGYDEPOSIT_CXX__

#include "SuperaSimEnergyDeposit.h"
#include "GenRandom.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaSimEnergyDepositProcessFactory __global_SuperaSimEnergyDepositProcessFactory__;

  SuperaSimEnergyDeposit::SuperaSimEnergyDeposit(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaSimEnergyDeposit::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _particle_label = cfg.get<std::string>("ParticleProducer");
    _output_label   = cfg.get<std::string>("OutCluster3DLabel");
    _store_dx = cfg.get<bool>("StoreLength");
    _store_dq = cfg.get<bool>("StoreCharge");
    _store_dp = cfg.get<bool>("StorePhoton");
    _store_dt = cfg.get<bool>("StoreDiffTime");
    _store_at = cfg.get<bool>("StoreAbsTime");
    _store_dedx = cfg.get<bool>("StoreDEDX");
 
  }

  void SuperaSimEnergyDeposit::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaSimEnergyDeposit::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    // Main cluster3d to be filled
    auto& event_cluster_v = mgr.get_data<larcv::EventClusterVoxel3D>(_output_label);
    auto const& meta = event_cluster_v.meta();

    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();

    std::vector<larcv::VoxelSet> cluster_de_v;
    std::vector<larcv::VoxelSet> cluster_dx_v;
    std::vector<larcv::VoxelSet> cluster_dq_v;
    std::vector<larcv::VoxelSet> cluster_dp_v;
    std::vector<larcv::VoxelSet> cluster_dt_v;
    std::vector<larcv::VoxelSet> cluster_at_v;
    std::vector<larcv::VoxelSet> cluster_dedx_v;

    // List particles to be stored
    static std::vector<int> part_idx_v(1e6,-1);
    std::fill(part_idx_v.begin(),part_idx_v.end(),-1);    
    auto const& mcp_v = mgr.get_data<larcv::EventParticle>(_particle_label).as_vector();
    LARCV_INFO() << "Processing larcv::EventParticle array: " << mcp_v.size() << std::endl;
    for(size_t idx=0; idx<mcp_v.size(); ++idx) {
      auto const& mcp = mcp_v[idx];
      auto const track_id = mcp.track_id();
      LARCV_INFO() << "Recording MCParticle track ID " << track_id << std::endl;
      if(part_idx_v.size() <= track_id) part_idx_v.resize(track_id+1,-1);
      part_idx_v[track_id] = idx;
    }

    // Reserve size
    cluster_de_v.resize(mcp_v.size()+1);
    if(_store_dq) cluster_dq_v.resize(mcp_v.size()+1);
    if(_store_dp) cluster_dp_v.resize(mcp_v.size()+1);
    if(_store_dt) cluster_dt_v.resize(mcp_v.size()+1);
    if(_store_at) cluster_at_v.resize(mcp_v.size()+1);
    if(_store_dx || _store_dedx) cluster_dx_v.resize(mcp_v.size()+1);
    // Register particle energy deposition coordinates
    auto const& sedep_v = LArData<supera::LArSimEnergyDeposit_t>();
    LARCV_INFO() << "Processing SimEnergyDeposit array: " << sedep_v.size() << std::endl;
    std::vector<int> store_idx_v;
    store_idx_v.reserve(sedep_v.size());
    for(size_t sedep_idx=0; sedep_idx<sedep_v.size(); ++sedep_idx) {
      auto const& sedep = sedep_v.at(sedep_idx);

      larcv::Point3D pt;
      VoxelID_t vox_id = meta.id(sedep.X(), sedep.Y(), sedep.Z());
      if(vox_id == larcv::kINVALID_VOXELID) {
	LARCV_DEBUG() << "Skipping sedep from track id " << sedep.TrackID() 
		      << " E=" << sedep.Energy()
		      << " pos=(" << sedep.X() << "," << sedep.Y() << "," << sedep.Z() << ")" << std::endl;
	continue;
      }
      LARCV_DEBUG() << "Recording sedep from track id " << sedep.TrackID() 
		    << " E=" << sedep.Energy() << std::endl;
      size_t cluster_idx = mcp_v.size();
      int track_id = sedep.TrackID();
      if(track_id < 0) track_id *= -1;
      if(track_id < (int)(part_idx_v.size()) && part_idx_v[track_id]>=0)
	cluster_idx = part_idx_v[track_id];

      float de = sedep.Energy();
      float dx = sedep.StepLength();
      float dq = sedep.NumElectrons();
      float dp = sedep.NumPhotons();
      float at = sedep.T();
      float dt = sedep.EndT() - sedep.StartT();
      
      cluster_de_v[cluster_idx].emplace(vox_id, de, true);

      if(_store_dq) cluster_dq_v[cluster_idx].emplace  (vox_id, dq, true);
      if(_store_dp) cluster_dp_v[cluster_idx].emplace  (vox_id, dp, true);
      if(_store_dt) cluster_dt_v[cluster_idx].emplace  (vox_id, dt, true);
      if(_store_dx || _store_dedx) cluster_dx_v[cluster_idx].emplace  (vox_id, dx, true);
      if(_store_at) {
	auto& cluster = cluster_at_v[cluster_idx];
	auto const& vox = cluster.find(vox_id);
	if(vox.id() == larcv::kINVALID_VOXELID)
	  cluster.emplace(vox_id, at, true);
	else if(at < vox.value())
	  cluster.emplace(vox_id,at,false);
      }
    }

    if(_store_dedx) {
      cluster_dedx_v.resize(mcp_v.size()+1);
      for(size_t cluster_idx=0; cluster_idx<cluster_de_v.size(); ++cluster_idx) {
	auto const& cluster_de = cluster_de_v[cluster_idx].as_vector();
	auto const& cluster_dx = cluster_dx_v[cluster_idx].as_vector();
	auto& cluster_dedx     = cluster_dedx_v[cluster_idx];
	for(size_t vox_idx=0; vox_idx < cluster_de.size(); ++vox_idx) {
	  auto const& vox_de = cluster_de[vox_idx];
	  auto const& vox_dx = cluster_dx[vox_idx];
	  float dedx = vox_de.value() / vox_dx.value();
	  //std::cout<< vox_de.value() << " " << vox_dx.value() << " " << dedx << std::endl;
	  cluster_dedx.emplace(vox_de.id(),dedx,true);
	}
      }
    }
    larcv::VoxelSetArray vsa_de; vsa_de.emplace(std::move(cluster_de_v));
    larcv::VoxelSetArray vsa_dx; vsa_dx.emplace(std::move(cluster_dx_v));
    larcv::VoxelSetArray vsa_dq; vsa_dq.emplace(std::move(cluster_dq_v));
    larcv::VoxelSetArray vsa_dp; vsa_dp.emplace(std::move(cluster_dp_v));
    larcv::VoxelSetArray vsa_dt; vsa_dt.emplace(std::move(cluster_dt_v));
    larcv::VoxelSetArray vsa_at; vsa_at.emplace(std::move(cluster_at_v));
    larcv::VoxelSetArray vsa_dedx; vsa_dedx.emplace(std::move(cluster_dedx_v));

    event_cluster_v.emplace(std::move(vsa_de),meta);
    if(_store_dx) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dx").emplace(std::move(vsa_dx),meta);
    if(_store_dq) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dq").emplace(std::move(vsa_dq),meta);
    if(_store_dp) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dp").emplace(std::move(vsa_dp),meta);
    if(_store_dt) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dt").emplace(std::move(vsa_dt),meta);
    if(_store_at) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_at").emplace(std::move(vsa_at),meta);
    if(_store_dedx) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dedx").emplace(std::move(vsa_dedx),meta);
    
    return true;
  }
      
  void SuperaSimEnergyDeposit::finalize()
  {}

}

#endif
