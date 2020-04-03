#ifndef __SuperaTrue2RecoVoxel3D_CXX__
#define __SuperaTrue2RecoVoxel3D_CXX__

#include <unordered_set>
#include <fstream>
#include <algorithm>
#include "SuperaTrue2RecoVoxel3D.h"
#include "RecoVoxel3D.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

template <typename T, typename M>
inline void dump_sim_channels(T const& sim_chnls, M const& meta, std::string fname)
{
  std::ofstream out(fname);
  out << "ch,time,track_id,n_e,energy,x,y,z,id\n";
  for (auto const& sim_ch: sim_chnls) {
    auto ch = sim_ch.Channel();

    for (auto const tick_ides : sim_ch.TDCIDEMap()) {
      double time = supera::TPCTDC2Tick(tick_ides.first);
	    for (auto const& edep : tick_ides.second) {
        out << ch << ','
          << time << ','
          << edep.trackID << ','
          << edep.numElectrons << ','
          << edep.energy << ','
          << edep.x << ','
          << edep.y << ','
          << edep.z << ','
          << meta.id(edep.x, edep.y, edep.z) << '\n';
      }
    }
  }
  out.close();
}

template <typename T>
inline void dump_reco_hits(T const& hits, std::string fname)
{
  std::ofstream out(fname);
  out << "ch,time,rms,amp,charge\n";

  for (auto const& hit : hits) {
    out << hit.Channel() << ","
      << hit.PeakTime() << ","
      << hit.RMS() << ","
      << hit.PeakAmplitude() << ","
      << hit.Integral() << "\n";
  }
  out.close();
}

template <typename T, typename F, typename M>
inline void dump_cluster3d(
    T const& space_pts,
    F const& finder,
    M const& meta,
    std::string fname)
{
  std::ofstream out(fname);
  out << "id,x,y,z,charge,ch1,t1,rms1,ch2,t2,rms2,ch3,t3,rms3\n";

  for (size_t i = 0; i < space_pts.size(); ++i) {
    auto const &pt = space_pts[i];
    auto *xyz = pt.XYZ();

    out
      << meta.id(xyz[0], xyz[1], xyz[2]) << ','
      << xyz[0] << ','
      << xyz[1] << ','
      << xyz[2] << ','
      << pt.ErrXYZ()[1];

    std::vector<art::Ptr<recob::Hit>> hits;
    finder.get(i, hits);

    for (auto const& hit : hits) {
      out << ','
        << hit->Channel() << ','
        << hit->PeakTime() << ','
        << hit->RMS();
    }
    out << "\n";
  }
  out.close();
}

template <typename T, typename M>
inline void dump_ghosts(T const& ghosts, M const& meta, std::string fname)
{
  std::ofstream out(fname);
  out << "id,x,y,z,is_ghost\n";

  for (const auto& [vox_id, is_ghost]: ghosts) {
    out
      << vox_id << ','
      << meta.pos_x(vox_id) << ','
      << meta.pos_y(vox_id) << ','
      << meta.pos_z(vox_id) << ','
      << is_ghost << '\n';
  }
  out.close();
}

template <typename T, typename M>
inline void dump_voxels(T const& voxels, M const& meta, std::string fname)
{
  std::ofstream out(fname);
  out << "id,x,y,z\n";

  for (const auto& vox_id : voxels) {
    out
      << vox_id << ','
      << meta.pos_x(vox_id) << ','
      << meta.pos_y(vox_id) << ','
      << meta.pos_z(vox_id) << '\n';
  }
  out.close();
}

template <typename T>
inline void dump_reco2true(T const& reco2true, std::string fname)
{
  std::ofstream out(fname);
  out << "reco_id,true_id,track_id\n";
  for (auto const& [reco_id, true_hits] : reco2true) {
    for (auto const& hit : true_hits) {
      out
        << reco_id << ','
        << hit.voxel_id << ','
        << hit.track_id << '\n';
    }
  }
  out.close();
}

template <typename T>
inline void dump_true2reco(T const& true2reco, std::string fname)
{
  std::ofstream out(fname);
  out << "true_id,reco_id\n";
  for (auto& [true_id, reco_pts] : true2reco) {
    for (auto reco_pt : reco_pts) {
      auto reco_id = reco_pt.get_id();
      out
        << true_id << ','
        << reco_id << '\n';
    }
  }
  out.close();
}

namespace larcv {

  static SuperaTrue2RecoVoxel3DProcessFactory __global_SuperaTrue2RecoVoxel3DProcessFactory__;

  SuperaTrue2RecoVoxel3D::SuperaTrue2RecoVoxel3D(const std::string name)
    : SuperaBase(name)
  {}

  void SuperaTrue2RecoVoxel3D::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _output_tensor3d = cfg.get<std::string>("OutputTensor3D","");
    _output_cluster3d = cfg.get<std::string>("OutputCluster3D","");
    _debug = cfg.get<bool> ("DebugMode",false);
    _hit_producer = cfg.get<std::string>("LArHitProducer","gaushit");
    _sps_producer = cfg.get<std::string>("LArSpacePointProducer","cluster3d");
    _use_true_pos = cfg.get<bool>("UseTruePosition",true);
    _twofold_matching = cfg.get<bool>("TwofoldMatching", false);
    _ref_meta3d_cluster3d = cfg.get<std::string>("Meta3DFromCluster3D","pcluster");
    _ref_meta3d_tensor3d = cfg.get<std::string>("Meta3DFromTensor3D","");
    _hit_threshold_ne = cfg.get<double>("HitThresholdNe", 0);
    _hit_window_ticks = cfg.get<double>("HitWindowTicks", 5.);
    _hit_peak_finding = cfg.get<bool>("HitPeakFinding", false);
    _dump_to_csv = cfg.get<bool>("DumpToCSV", false);
    _post_averaging = cfg.get<bool>("PostAveraging", false);
    _post_averaging_threshold = cfg.get<double>("PostAveragingThreshold_cm", 0.3);
  }


  void SuperaTrue2RecoVoxel3D::initialize()
  {
    SuperaBase::initialize();

    // dump channel map
    if (_dump_to_csv) {
      art::ServiceHandle<geo::Geometry const> geom;

      std::ofstream out("channel_map.csv");
      out << "ch,cryo,tpc,plane,wire\n";

      for (size_t ch = 0; ch < geom->Nchannels(); ++ch) {

        auto const& ids = geom->ChannelToWire(ch);
        for (auto const& id : ids) {
          out << ch << ','
            << id.Cryostat << ','
            << id.TPC << ','
            << id.Plane << ','
            << id.Wire << '\n';
        }
      }
      out.close();
    }
  }

  larcv::Voxel3DMeta SuperaTrue2RecoVoxel3D::get_meta3d(IOManager& mgr) const {

    larcv::Voxel3DMeta meta3d;
    if(!_ref_meta3d_cluster3d.empty()) {
      auto const& ev_cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(_ref_meta3d_cluster3d);
      meta3d = ev_cluster3d.meta();
    }
    else if(!_ref_meta3d_tensor3d.empty()) {
      auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_ref_meta3d_tensor3d);
      meta3d = ev_tensor3d.meta();
    }
    else
      LARCV_CRITICAL() << "ref_meta3d_cluster3d nor ref_meta3d_tensor3d set!" << std::endl;
    return meta3d;
  }

  bool SuperaTrue2RecoVoxel3D::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);
    //
    // Step 0. get 3D meta
    // Step 1. create a list of true hits
    //   - store an array of hits per channel, each hit contains time and 3D voxel ID)
    // Step 2. Create a map of reco 3D voxel <=> true 3D voxel (per plane, as a representation of matched hits)
    // Step 3. Identify valid reco3D voxels
    //   - Find an overlapping true 3D voxel ID across planes per 3D reco voxel ID
    //

    //
    // Step 0. ... get meta3d
    //
		auto event_id = mgr.event_id().event();
    LARCV_INFO() << "Retrieving 3D meta..." << std::endl;
    auto meta3d = get_meta3d(mgr);

		std::unordered_set<VoxelID_t> true_voxel_ids;

    //
    // Step 1. ... create a list of true hits
    //

    // Get geometry info handler
    auto geop = lar::providerFrom<geo::Geometry>();

    // Create a hit list container
    // true_hit_vv[ch][i_hit]
    std::vector<std::vector<TrueHit_t> > true_hit_vv;
    true_hit_vv.resize(geop->Nchannels());

    // Fill a hit list container
    for(auto const& sch : LArData<supera::LArSimCh_t>()){
      // Get a unique readout channel number
      auto ch = sch.Channel();

      // Loop over hits and store
      for (auto const tick_ides : sch.TDCIDEMap()) {
        double x_pos = (supera::TPCTDC2Tick(tick_ides.first) * supera::TPCTickPeriod() - supera::TriggerOffsetTPC()) * supera::DriftVelocity();
        TrueHit_t hit;
        hit.time = supera::TPCTDC2Tick(tick_ides.first);

        for (auto const& edep : tick_ides.second) {
          if (edep.numElectrons < _hit_threshold_ne) continue;

          if(_use_true_pos) x_pos = edep.x;

          auto vox_id = meta3d.id(x_pos, edep.y, edep.z);
      	  if(vox_id == larcv::kINVALID_VOXELID) continue;
          true_voxel_ids.insert(vox_id);

          hit.track_voxel_ids.emplace_back(vox_id, edep.trackID);
          hit.n_electrons.push_back(edep.numElectrons);
        }

        // append hit to true_hit_vv[ch]
        if(hit.track_voxel_ids.size())
          true_hit_vv[ch].push_back(std::move(hit));
      }
    }

    LARCV_INFO() << "Created a list of true hits: " << true_hit_vv.size() << " channels" << std::endl;
    if(_debug) {
      size_t num_hits = 0;
      for(auto const& true_hit_v : true_hit_vv) num_hits += true_hit_v.size();
      LARCV_INFO() << " ... corresponds to " << num_hits << " true hits!" << std::endl;
    }

    // ---------------------------------------------------
    // Loop over 3d space point
    // For each 3d point, find the associated reco 2d hits
    // Match reco 2d hit -> sim hit
    // Check whether all sim hits are originated from the
    // same ionization position (in voxels)
    // ---------------------------------------------------

    auto const *ev = GetEvent();
    auto space_pts = ev->getValidHandle<std::vector<recob::SpacePoint>>(_sps_producer);

    // TODO(2020-03-20 kvtsang)  No space point, warnning?
    if (!space_pts.isValid()) return true;
    art::FindManyP<recob::Hit> hit_finder(space_pts, *ev, _sps_producer);

    // reco2true : reco_voxel_id -> [(true_voxel_id, track_id)]
    std::map<VoxelID_t, std::set<TrackVoxel_t>> reco2true;

    for (size_t i = 0; i < space_pts->size(); ++i) {
      auto const &pt = space_pts->at(i);
      auto *xyz = pt.XYZ();

      auto reco_voxel_id = meta3d.id(xyz[0], xyz[1], xyz[2]);
      if(reco_voxel_id == larcv::kINVALID_VOXELID) continue;

      std::vector<art::Ptr<recob::Hit>> hits;
      hit_finder.get(i, hits);

      // matching: gaushit -> simhit -> true 3d voxel
      // 3 voxel set per 1 space point
      std::vector<std::set<TrackVoxel_t>> matched_voxels;

      for (auto const& hit_ptr : hits) {
        auto const& hit = *hit_ptr;

        auto ch = hit.Channel();
        auto t = hit.PeakTime();

        double t_start = t - _hit_window_ticks;
        double t_end = t + _hit_window_ticks;

        std::set<TrackVoxel_t> track_voxel_ids;
        if (_hit_peak_finding)
          find_hit_peaks(true_hit_vv[ch], t_start, t_end, track_voxel_ids);
        else
          find_hits(true_hit_vv[ch], t_start, t_end, track_voxel_ids);

        matched_voxels.push_back(std::move(track_voxel_ids));
      }

      // finding overlaps
      std::set<TrackVoxel_t> overlaps;

      if (_twofold_matching) {
        // Optional: require matches from two different plane
        for (size_t p1 = 0; p1 < matched_voxels.size(); ++p1) {
          auto const& v1 = matched_voxels[p1];
          for (size_t p2 = p1 + 1; p2 < matched_voxels.size(); ++p2) {
            auto const& v2 = matched_voxels[p2];

            //overlaps between planes p1 and p2
            std::set_intersection(
                v1.begin(), v1.end(),
                v2.begin(), v2.end(),
                std::inserter(overlaps, overlaps.end()));
          }
        }
      }
      else {
        // non-ghost = if all hits from different planes come from the same true voxel
        if (matched_voxels.size() > 0) {
          auto& v = matched_voxels[0];
          overlaps.insert(v.begin(), v.end());
        }

        for (size_t plane = 1; plane < matched_voxels.size(); ++plane) {
            auto const& v = matched_voxels[plane];

            // temporary storage
            std::set<TrackVoxel_t> overlaps_;

            std::set_intersection(
                overlaps.begin(), overlaps.end(),
                v.begin(), v.end(),
                std::inserter(overlaps_, overlaps_.end()));

            overlaps = std::move(overlaps_);
        }
      }

      auto itr = reco2true.find(reco_voxel_id);
      if (itr == reco2true.end())
        reco2true[reco_voxel_id] = std::move(overlaps);
      else
        itr->second.insert(overlaps.begin(), overlaps.end());

    } // end looping reco pts

    // FIXME(2020-03-20 kvtsang) output ghost point labels to VoxelSet
    std::map<VoxelID_t, bool> ghosts;
    for (auto const& key_value : reco2true) ghosts.emplace(key_value.first, true);

    // true2reco : true_voxel_id -> [reco_voxel_id]
    // inverting reco2true for non-ghost points, track_id contracted
    std::map<VoxelID_t, std::unordered_set<RecoVoxel3D>> true2reco;
    auto insert_true2reco = [&](VoxelID_t true_id, VoxelID_t reco_id) {
      RecoVoxel3D reco_pt(reco_id);

      auto itr = true2reco.find(true_id);
      if (itr == true2reco.end())
        true2reco.emplace(true_id, std::unordered_set<RecoVoxel3D>({reco_pt}));
      else
        itr->second.insert(reco_pt);
    };

    if (_post_averaging)
      set_ghost_with_averaging(reco2true, meta3d, ghosts, insert_true2reco);
    else
      set_ghost(reco2true, ghosts, insert_true2reco);

    LARCV_INFO() << true2reco.size() << " true points mapped to " << reco2true.size() << " reco points" << std::endl;
    // store true energy in VoxelSet
    if(!_output_tensor3d.empty() || !_output_cluster3d.empty()) {
      EventSparseTensor3D* event_tensor3d = nullptr;
      EventClusterVoxel3D* event_cluster3d = nullptr;
      if(!_output_tensor3d.empty()) {
        event_tensor3d = (larcv::EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_tensor3d));
        event_tensor3d->reserve(true2reco.size());
      }
      if(!_output_cluster3d.empty()) {
        event_cluster3d = (larcv::EventClusterVoxel3D*)(mgr.get_data("cluster3d",_output_cluster3d));
        event_cluster3d->resize(true2reco.size());
      }
      
      size_t cluster_ctr=0;
      for(auto const& keyval : true2reco) {
        if(event_cluster3d && keyval.second.empty()) continue;
        if(event_tensor3d) event_tensor3d->emplace(keyval.first, 0., true);
        auto& vs = event_cluster3d->writeable_voxel_set(cluster_ctr);
        vs.reserve(keyval.second.size());
        for(auto const& pt: keyval.second) vs.emplace(pt.get_id(), 0., true);
        ++cluster_ctr;
      }
    }

    // store corresponding reco points in VoxelSetArray (outer index == true VoxelSet index)

    // --------------------------------
    //  ____  _____ ____  _   _  ____
    //  |  _ \| ____| __ )| | | |/ ___|
    //  | | | |  _| |  _ \| | | | |  _
    //  | |_| | |___| |_) | |_| | |_| |
    //  |____/|_____|____/ \___/ \____|
    // --------------------------------

    if (!_dump_to_csv) return true;

    auto save_to = [&](std::string prefix) {
      return prefix + "_" + std::to_string(event_id) + ".csv";
    };

    // SimChannel
    dump_sim_channels(LArData<supera::LArSimCh_t>(), meta3d, save_to("simch"));

    // cluster3d
    dump_cluster3d(*space_pts, hit_finder, meta3d, save_to("reco3d"));

    // ghost label
    dump_ghosts(ghosts, meta3d, save_to("ghosts"));

    // true label (from SimChannels)
    dump_voxels(true_voxel_ids, meta3d, save_to("true3d"));

    // gaushit
    auto gaus_hits = ev->getValidHandle<std::vector<recob::Hit>>("gaushit");
    if (gaus_hits.isValid()) {
      dump_reco_hits(*gaus_hits, save_to("gaushit"));
    }

    // reco2true
    dump_reco2true(reco2true, save_to("reco2true"));

    // true2reco_nonghost (after averaging)
    dump_true2reco(true2reco, save_to("true2reco"));
    return true;
  }


  void SuperaTrue2RecoVoxel3D::set_ghost_with_averaging(
      const std::map<VoxelID_t, std::set<TrackVoxel_t>>& reco2true,
      const larcv::Voxel3DMeta& meta3d,
      std::map<VoxelID_t, bool>& ghosts,
      std::function<void(VoxelID_t, VoxelID_t)> const& insert_true2reco)
  {
    double threshold2 = pow(_post_averaging_threshold, 2);

    // build inverse map (true_voxel_id, track_id) -> [reco_voxel_id]
    // use later for averaging
    std::map<TrackVoxel_t, std::unordered_set<VoxelID_t>> true2reco;
    for (auto const& [reco_id, true_hits] : reco2true) {
      for (auto const& hit : true_hits) {
          auto itr = true2reco.find(hit);
          if (itr == true2reco.end())
            true2reco.emplace(hit, std::unordered_set<VoxelID_t>({reco_id}));
          else
            itr->second.insert(reco_id);
      } // loop true hits
    } // loop reco2true

    // set ghost label, and
    for (auto const& [true_hit, reco_voxels] : true2reco) {
      size_t n = reco_voxels.size();

      // 1-to-1 match
      // mark as non-ghost
      if (n == 1) {
        auto reco_id = *reco_voxels.begin();
        ghosts[reco_id] = false;
        insert_true2reco(true_hit.voxel_id, reco_id);
        continue;
      }

      // 1-to-many match
      // calcuate the mean reco. position (per true voxel + track_id)
      // mark a subset as non-ghost near mean
      std::vector<VoxelID_t> ids;
      std::vector<double> x, y, z;

      // convert reco voxel ids to xyz
      for (auto id: reco_voxels) {
        ids.push_back(id);
        x.push_back(meta3d.pos_x(id));
        y.push_back(meta3d.pos_y(id));
        z.push_back(meta3d.pos_z(id));
      }

      // mean reco position for a given true pt
      double x0 = std::accumulate(x.cbegin(), x.cend(), 0.) / n;
      double y0 = std::accumulate(y.cbegin(), y.cend(), 0.) / n;
      double z0 = std::accumulate(z.cbegin(), z.cend(), 0.) / n;

      // keep reco pts with a threshold arround reco mean position
      for (size_t i = 0; i< n; ++i) {
        auto reco_id = ids[i];

        double dx = x[i] - x0;
        double dy = y[i] - y0;
        double dz = z[i] - z0;
        double dist2 = dx*dx + dy*dy + dz*dz;

        if (dist2 < threshold2) {
          ghosts[reco_id] = false;
          insert_true2reco(true_hit.voxel_id, reco_id);
        }
      } // average over reco pts
    } // loop true2reco
  }

  void SuperaTrue2RecoVoxel3D::set_ghost(
        const std::map<VoxelID_t, std::set<TrackVoxel_t>>& reco2true,
        std::map<VoxelID_t, bool>& ghosts,
        std::function<void(VoxelID_t, VoxelID_t)> const& insert_true2reco)
  {
    for (auto const& [reco_id, true_hits] : reco2true) {
      if (true_hits.size() == 0) continue; // no overlapping hit
      ghosts[reco_id] = false;
      for (auto const& hit : true_hits) insert_true2reco(hit.voxel_id, reco_id);
    }
  }

  void SuperaTrue2RecoVoxel3D::find_hit_peaks(const std::vector<TrueHit_t>& hits,
      double t_start, double t_end, std::set<TrackVoxel_t>& track_voxel_ids)
  {
    std::map<int, size_t> track_idx;
    std::vector<double> n_electrons;
    std::vector<VoxelID_t> voxel_ids;

    auto update = [&](int track_id, VoxelID_t voxel_id, double ne) {
      auto itr = track_idx.find(track_id);
      if (itr == track_idx.end()) {
        size_t idx = track_idx.size();
        track_idx.emplace(track_id, idx);
        n_electrons.push_back(ne);
        voxel_ids.push_back(track_id);
      }
      else {
        size_t idx = itr->second;
        if (ne > n_electrons[idx]) {
          n_electrons[idx] = ne;
          voxel_ids[idx] = voxel_id;
        }
      }
    };

    for (auto const& hit: hits) {
      if (t_start <= hit.time && hit.time <= t_end) {
        for (size_t i =0; i < hit.track_voxel_ids.size(); ++i) {
          auto const& true_info = hit.track_voxel_ids[i];
          update(true_info.track_id, true_info.voxel_id, hit.n_electrons[i]);

        }// loop (track_id, voxel_id)
      } // hit time selection
    } // loop hits

    // insert (voxel_id, track_id) from peaks
    for (auto [track_id, idx] : track_idx) {
      track_voxel_ids.emplace(voxel_ids[idx], track_id);
    }
  }

  void SuperaTrue2RecoVoxel3D::find_hits(const std::vector<TrueHit_t>& hits,
      double t_start, double t_end, std::set<TrackVoxel_t>& track_voxel_ids)
  {
    for (auto const& hit: hits) {
      // TODO(2020-03-02 kvtsang) binary search?
      // Hit times are sorted.
      // In principle could use binary search if need better performance
      if (t_start <= hit.time && hit.time <= t_end)
        track_voxel_ids.insert(
            hit.track_voxel_ids.begin(),
            hit.track_voxel_ids.end());
    }
  }

  void SuperaTrue2RecoVoxel3D::finalize()
  {}
}

#endif
