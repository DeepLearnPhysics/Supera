#ifndef __SUPERASIMENERGYDEPOSIT_CXX__
#define __SUPERASIMENERGYDEPOSIT_CXX__

#include "SuperaSimEnergyDeposit.h"
#include "GenRandom.h"

namespace larcv {

  static SuperaSimEnergyDepositProcessFactory __global_SuperaSimEnergyDepositProcessFactory__;

  SuperaSimEnergyDeposit::SuperaSimEnergyDeposit(const std::string name) : SuperaBase(name) {}
  void SuperaSimEnergyDeposit::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _particle_label = cfg.get<std::string>("ParticleProducer");
    _output_label = cfg.get<std::string>("OutCluster3DLabel");
    _store_dx = cfg.get<bool>("StoreLength");
    _store_dq = cfg.get<bool>("StoreCharge");
    _store_dp = cfg.get<bool>("StorePhoton");
    _store_dt = cfg.get<bool>("StoreDiffTime");
    _store_at = cfg.get<bool>("StoreAbsTime");
    _store_dedx = cfg.get<bool>("StoreDEDX");
    _use_lite = cfg.get<bool>("UseSEDLite", false);

    auto cryostat_v = cfg.get<std::vector<unsigned short>>("CryostatList");
    auto tpc_v = cfg.get<std::vector<unsigned short>>("TPCList");
    assert(cryostat_v.size() == tpc_v.size());
    larcv::Point3D min_pt(1.e9, 1.e9, 1.e9);
    larcv::Point3D max_pt(-1.e9, -1.e9, -1.e9);
    auto geop = lar::providerFrom<geo::Geometry>();
    for (size_t idx = 0; idx < cryostat_v.size(); ++idx) {
      auto const& c = cryostat_v[idx];
      auto const& t = tpc_v[idx];
      auto const& cryostat = geop->Cryostat(geo::CryostatID(c));
      if (!cryostat.HasTPC(t)) {
        LARCV_CRITICAL() << "Invalid TPCList: cryostat " << c << " does not contain tpc " << t
                         << std::endl;
        throw larbys();
      }
      auto const& tpcabox = cryostat.TPC(t).ActiveBoundingBox();
      if (min_pt.x > tpcabox.MinX()) min_pt.x = tpcabox.MinX();
      if (min_pt.y > tpcabox.MinY()) min_pt.y = tpcabox.MinY();
      if (min_pt.z > tpcabox.MinZ()) min_pt.z = tpcabox.MinZ();
      if (max_pt.x < tpcabox.MaxX()) max_pt.x = tpcabox.MaxX();
      if (max_pt.y < tpcabox.MaxY()) max_pt.y = tpcabox.MaxY();
      if (max_pt.z < tpcabox.MaxZ()) max_pt.z = tpcabox.MaxZ();
    }
    _world_bounds.update(min_pt, max_pt);
  }

  void SuperaSimEnergyDeposit::initialize()
  {
    SuperaBase::initialize();
  }

  void SuperaSimEnergyDeposit::fill_sedep_v(const std::vector<larcv::Particle> mcp_v,
                                            std::vector<int>& part_idx_v,
                                            larcv::Voxel3DMeta meta,
                                            larcv::EventClusterVoxel3D& event_de_v)
  {

    auto const& sedep_v = LArData<supera::LArSimEnergyDepositLite_t>();
    LARCV_INFO() << "Processing SimEnergyDeposit array: " << sedep_v.size() << std::endl;
    for (size_t sedep_idx = 0; sedep_idx < sedep_v.size(); ++sedep_idx) {
      auto const& sedep = sedep_v.at(sedep_idx);

      larcv::Point3D pt;
      VoxelID_t vox_id = meta.id(sedep.X(), sedep.Y(), sedep.Z());
      int track_id = std::abs(sedep.TrackID());
      if (vox_id == larcv::kINVALID_VOXELID ||
          !_world_bounds.contains(sedep.X(), sedep.Y(), sedep.Z())) {
        LARCV_DEBUG() << "Skipping sedep from track id " << track_id << " E=" << sedep.Energy()
                      << " pos=(" << sedep.X() << "," << sedep.Y() << "," << sedep.Z() << ")"
                      << std::endl;
        continue;
      }
      LARCV_DEBUG() << "Recording sedep from track id " << track_id << " E=" << sedep.Energy()
                    << std::endl;
      size_t cluster_idx = mcp_v.size();
      if (track_id < (int)(part_idx_v.size()) && part_idx_v.at(track_id) >= 0)
        cluster_idx = part_idx_v.at(track_id);
      assert(cluster_idx == event_de_v.as_vector().size());
      LARCV_DEBUG() << "Cluster index: " << cluster_idx << " / " << event_de_v.as_vector().size()
                    << std::endl;
      float de = sedep.Energy();
      //float at = sedep.T();

      larcv::Voxel v(vox_id, de);
      event_de_v.writeable_voxel_set(cluster_idx).add(v);
    } // end for
  }

  void SuperaSimEnergyDeposit::fill_sedep_v(const std::vector<larcv::Particle> mcp_v,
                                            std::vector<int>& part_idx_v,
                                            larcv::Voxel3DMeta meta,
                                            larcv::EventClusterVoxel3D& event_de_v,
                                            larcv::EventClusterVoxel3D* event_dq_v,
                                            larcv::EventClusterVoxel3D* event_dp_v,
                                            larcv::EventClusterVoxel3D* event_dx_v,
                                            larcv::EventClusterVoxel3D* event_dt_v,
                                            larcv::EventClusterVoxel3D* event_at_v,
                                            larcv::EventClusterVoxel3D* event_dedx_v)
  {
    auto const& sedep_v = LArData<supera::LArSimEnergyDeposit_t>();
    LARCV_INFO() << "Processing SimEnergyDeposit array: " << sedep_v.size() << std::endl;
    for (size_t sedep_idx = 0; sedep_idx < sedep_v.size(); ++sedep_idx) {
      auto const& sedep = sedep_v.at(sedep_idx);

      larcv::Point3D pt;
      VoxelID_t vox_id = meta.id(sedep.X(), sedep.Y(), sedep.Z());
      int track_id = std::abs(sedep.TrackID());
      if (vox_id == larcv::kINVALID_VOXELID ||
          !_world_bounds.contains(sedep.X(), sedep.Y(), sedep.Z())) {
        LARCV_DEBUG() << "Skipping sedep from track id " << track_id << " E=" << sedep.Energy()
                      << " pos=(" << sedep.X() << "," << sedep.Y() << "," << sedep.Z() << ")"
                      << std::endl;
        continue;
      }
      LARCV_DEBUG() << "Recording sedep from track id " << track_id << " E=" << sedep.Energy()
                    << std::endl;
      size_t cluster_idx = mcp_v.size();
      if (track_id < (int)(part_idx_v.size()) && part_idx_v.at(track_id) >= 0)
        cluster_idx = part_idx_v.at(track_id);
      assert(cluster_idx == event_de_v.as_vector().size());
      LARCV_DEBUG() << "Cluster index: " << cluster_idx << " / " << event_de_v.as_vector().size()
                    << std::endl;
      float de = sedep.Energy();
      float at = sedep.T();

      larcv::Voxel v(vox_id, de);
      event_de_v.writeable_voxel_set(cluster_idx).add(v);

      // What follows is only to store extra information
      // (afaik not used currently)

      if (event_dq_v != nullptr) {
        assert(event_dp_v != nullptr);
        assert(event_dx_v != nullptr);
        assert(event_dt_v != nullptr);
        assert(event_at_v != nullptr);
        assert(event_dedx_v != nullptr);

        assert(cluster_idx == event_dq_v->as_vector().size());
        assert(cluster_idx == event_dp_v->as_vector().size());
        assert(cluster_idx == event_dx_v->as_vector().size());
        assert(cluster_idx == event_dt_v->as_vector().size());
        assert(cluster_idx == event_at_v->as_vector().size());
        assert(cluster_idx == event_dedx_v->as_vector().size());
        if (_store_dq) {
          v.set(vox_id, sedep.NumElectrons());
          event_dq_v->writeable_voxel_set(cluster_idx).add(v);
        }
        if (_store_dp) {
          v.set(vox_id, sedep.NumPhotons());
          event_dp_v->writeable_voxel_set(cluster_idx).add(v);
        }
        if (_store_dt) {
          v.set(vox_id, sedep.EndT() - sedep.StartT());
          event_dt_v->writeable_voxel_set(cluster_idx).add(v);
        }
        if (_store_dx || _store_dedx) {
          v.set(vox_id, sedep.StepLength());
          event_dx_v->writeable_voxel_set(cluster_idx).add(v);
        }
        if (_store_at) {
          auto& cluster = event_at_v->writeable_voxel_set(cluster_idx);
          auto const& vox = cluster.find(vox_id);
          if (vox.id() == larcv::kINVALID_VOXELID)
            cluster.emplace(vox_id, at, true);
          else if (at < vox.value())
            cluster.emplace(vox_id, at, false);
        }
      }
    } //end for
  }   // end fill_sedep_v

  bool SuperaSimEnergyDeposit::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    // Main cluster3d to be filled
    auto& event_de_v = mgr.get_data<larcv::EventClusterVoxel3D>(_output_label);
    larcv::Voxel3DMeta meta = event_de_v.meta();

    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();
    /*
    std::vector<larcv::VoxelSet> cluster_de_v;
    std::vector<larcv::VoxelSet> cluster_dx_v;
    std::vector<larcv::VoxelSet> cluster_dq_v;
    std::vector<larcv::VoxelSet> cluster_dp_v;
    std::vector<larcv::VoxelSet> cluster_dt_v;
    std::vector<larcv::VoxelSet> cluster_at_v;
    std::vector<larcv::VoxelSet> cluster_dedx_v;
    */
    // List particles to be stored
    static std::vector<int> part_idx_v(1e6, -1);
    std::fill(part_idx_v.begin(), part_idx_v.end(), -1);
    auto const& mcp_v = mgr.get_data<larcv::EventParticle>(_particle_label).as_vector();
    LARCV_INFO() << "Processing larcv::EventParticle array: " << mcp_v.size() << std::endl;
    auto const& mcshower_v = LArData<supera::LArMCShower_t>();
    for (size_t idx = 0; idx < mcp_v.size(); ++idx) {
      auto const& mcp = mcp_v[idx];
      auto const track_id = mcp.track_id();
      LARCV_INFO() << "Recording MCParticle track ID " << track_id << " PDG " << mcp.pdg_code()
                   << std::endl;
      if (part_idx_v.size() <= track_id) part_idx_v.resize(track_id + 1, -1);
      part_idx_v.at(track_id) = idx;
      // If this is shower, look for a corresponding shower
      for (auto const& mcshower : mcshower_v) {
        if (mcshower.TrackID() != track_id) continue;
        for (auto const& daughter_id : mcshower.DaughterTrackID()) {
          if (part_idx_v.size() <= daughter_id) part_idx_v.resize(daughter_id + 1, -1);
          part_idx_v.at(daughter_id) = idx;
        }
        break;
      }
    }

    // Reserve size
    event_de_v.resize(mcp_v.size() + 1);
    larcv::EventClusterVoxel3D* event_dx_v = nullptr;
    larcv::EventClusterVoxel3D* event_dq_v = nullptr;
    larcv::EventClusterVoxel3D* event_dp_v = nullptr;
    larcv::EventClusterVoxel3D* event_dt_v = nullptr;
    larcv::EventClusterVoxel3D* event_at_v = nullptr;
    larcv::EventClusterVoxel3D* event_dedx_v = nullptr;

    if (_store_dx || _store_dedx) {
      event_dx_v = (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d", _output_label + "_dx"));
      event_dx_v->resize(mcp_v.size() + 1);
      event_dx_v->meta(meta);
    }
    if (_store_dq) {
      event_dq_v = (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d", _output_label + "_dq"));
      event_dq_v->resize(mcp_v.size() + 1);
      event_dq_v->meta(meta);
    }
    if (_store_dp) {
      event_dp_v = (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d", _output_label + "_dp"));
      event_dp_v->resize(mcp_v.size() + 1);
      event_dp_v->meta(meta);
    }
    if (_store_dt) {
      event_dt_v = (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d", _output_label + "_dt"));
      event_dt_v->resize(mcp_v.size() + 1);
      event_dt_v->meta(meta);
    }
    if (_store_at) {
      event_at_v = (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d", _output_label + "_at"));
      event_at_v->resize(mcp_v.size() + 1);
      event_at_v->meta(meta);
    }
    if (_store_dedx) {
      event_dedx_v =
        (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d", _output_label + "_dedx"));
      event_dedx_v->resize(mcp_v.size() + 1);
      event_dedx_v->meta(meta);
    }
    if (event_de_v.as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_de_v size mismatch with mcp_v: " << event_de_v.as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dx_v && event_dx_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dx_v size mismatch with mcp_v: " << event_dx_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dq_v && event_dq_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dq_v size mismatch with mcp_v: " << event_dq_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dp_v && event_dp_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dp_v size mismatch with mcp_v: " << event_dp_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dt_v && event_dt_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dt_v size mismatch with mcp_v: " << event_dt_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_at_v && event_at_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_at_v size mismatch with mcp_v: " << event_at_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dedx_v && event_dedx_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dedx_v size mismatch with mcp_v: " << event_dedx_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    // Register particle energy deposition coordinates

    if (_use_lite) { fill_sedep_v(mcp_v, part_idx_v, meta, event_de_v); }
    else {
      fill_sedep_v(mcp_v,
                   part_idx_v,
                   meta,
                   event_de_v,
                   event_dq_v,
                   event_dp_v,
                   event_dx_v,
                   event_dt_v,
                   event_at_v,
                   event_dedx_v);
    }

    if (_store_dedx) {
      LARCV_INFO() << "Computing dE/dX for clusters: dx " << event_dx_v->as_vector().size()
                   << " de " << event_de_v.as_vector().size() << std::endl;
      /*
      for(size_t cluster_idx=0; cluster_idx<cluster_de_v.size(); ++cluster_idx) {
	LARCV_INFO() << "Processing cluster " << cluster_idx << std::endl;
	//auto const& cluster_de = event_de_v.voxel_set(cluster_idx).as_vector();
	//auto const& cluster_dx = event_dx_v->voxel_set(cluster_idx).as_vector();
	auto const& cluster_de = event_de_v.as_vector().at(cluster_idx).as_vector();
	auto const& cluster_dx = event_dx_v->as_vector().at(cluster_idx).as_vector();
	auto& cluster_dedx     = event_dedx_v->writeable_voxel_set(cluster_idx);
	assert(cluster_de.size() == cluster_dx.size());
	for(size_t vox_idx=0; vox_idx < cluster_de.size(); ++vox_idx) {
	  auto const& vox_de = cluster_de[vox_idx];
	  auto const& vox_dx = cluster_dx[vox_idx];
	  float dedx = vox_de.value() / vox_dx.value();
	  cluster_dedx.emplace(vox_de.id(),dedx,true);
	}
      }
      */
    }

    // report
    for (size_t part_idx = 0; part_idx < mcp_v.size(); ++part_idx) {
      LARCV_INFO() << "Track ID " << mcp_v[part_idx].track_id() << " PDG "
                   << mcp_v[part_idx].pdg_code() << " Voxel Count "
                   << event_de_v.as_vector().at(part_idx).size() << std::endl;
    }

    // assert
    if (event_de_v.as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_de_v size mismatch with mcp_v: " << event_de_v.as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dx_v && event_dx_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dx_v size mismatch with mcp_v: " << event_dx_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dq_v && event_dq_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dq_v size mismatch with mcp_v: " << event_dq_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dp_v && event_dp_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dp_v size mismatch with mcp_v: " << event_dp_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dt_v && event_dt_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dt_v size mismatch with mcp_v: " << event_dt_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_at_v && event_at_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_at_v size mismatch with mcp_v: " << event_at_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    if (event_dedx_v && event_dedx_v->as_vector().size() != (mcp_v.size() + 1)) {
      LARCV_ERROR() << "event_dedx_v size mismatch with mcp_v: " << event_dedx_v->as_vector().size()
                    << " vs. " << mcp_v.size();
      throw std::exception();
    }
    /*
    larcv::VoxelSetArray vsa_de;   vsa_de.emplace(std::move(cluster_de_v));
    larcv::VoxelSetArray vsa_dx;   vsa_dx.emplace(std::move(cluster_dx_v));
    larcv::VoxelSetArray vsa_dq;   vsa_dq.emplace(std::move(cluster_dq_v));
    larcv::VoxelSetArray vsa_dp;   vsa_dp.emplace(std::move(cluster_dp_v));
    larcv::VoxelSetArray vsa_dt;   vsa_dt.emplace(std::move(cluster_dt_v));
    larcv::VoxelSetArray vsa_at;   vsa_at.emplace(std::move(cluster_at_v));
    larcv::VoxelSetArray vsa_dedx; vsa_dedx.emplace(std::move(cluster_dedx_v));

    event_cluster_v.set(vsa_de,meta);
    if(_store_dx)   mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dx"  ).set(vsa_dx,meta);
    if(_store_dq)   mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dq"  ).set(vsa_dq,meta);
    if(_store_dp)   mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dp"  ).set(vsa_dp,meta);
    if(_store_dt)   mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dt"  ).set(vsa_dt,meta);
    if(_store_at)   mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_at"  ).set(vsa_at,meta);
    if(_store_dedx) mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + "_dedx").set(vsa_dedx,meta);
    */
    return true;
  }

  void SuperaSimEnergyDeposit::finalize() {}

}

#endif
