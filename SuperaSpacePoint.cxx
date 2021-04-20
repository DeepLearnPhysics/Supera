#ifndef __SUPERASPACEPOINT_CXX__
#define __SUPERASPACEPOINT_CXX__

#include "SuperaSpacePoint.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

#include "canvas/Persistency/Common/FindManyP.h"

namespace larcv {

  static SuperaSpacePointProcessFactory __global_SuperaSpacePointProcessFactory__;

  SuperaSpacePoint::SuperaSpacePoint(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaSpacePoint::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _producer_labels    = cfg.get<std::vector<std::string>>("SpacePointProducers", {});
    _output_label       = cfg.get<std::string>("OutputLabel");
    _max_debug_dropping = cfg.get<size_t>("MaxDebugForDropping", 0);
    _n_planes           = cfg.get<size_t>("NumOfPlanes", 3);

    auto to_drop = cfg.get<std::vector<std::string>>("DropOutput", {});
    _drop_output.insert(to_drop.cbegin(), to_drop.cend());
    _store_wire_info = _drop_output.count("hit_*") == 0;

    _reco_charge_range = cfg.get<std::vector<double>>("RecoChargeRange", {0,9e99}); 
    assert(_reco_charge_range.size() == 2);

    // for backward compatibility
    auto prod_label = cfg.get<std::string>("SpacePointProducer", "");
    if (!prod_label.empty()) _producer_labels.push_back(prod_label);

    if (_producer_labels.size() == 0) {
      LARCV_ERROR() << "No SpacePointProducer(s)\n";
    }

    for (auto const& label : _producer_labels)
      Request(supera::LArDataType_t::kLArSpacePoint_t, label);
  }

  void SuperaSpacePoint::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaSpacePoint::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    // Get meta (presumbly produced by BBox)
    auto& main_tensor = mgr.get_data<larcv::EventSparseTensor3D>(_output_label);
    auto const& meta = main_tensor.meta();
    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();

    /* TODO(kvtsang) implement number of clusters
     * Now consider whole event as a single cluster
     */

    larcv::VoxelSet v_occupancy;
    larcv::VoxelSet v_charge;
    larcv::VoxelSet v_charge_asym;
    larcv::VoxelSet v_chi2;

    std::vector<larcv::VoxelSet> v_hit_charge(_n_planes);
    std::vector<larcv::VoxelSet> v_hit_amp   (_n_planes);
    std::vector<larcv::VoxelSet> v_hit_time  (_n_planes);
    std::vector<larcv::VoxelSet> v_hit_rms   (_n_planes);

    /* FIXME(kvtsang) To be removed?
     * Find associated hits
     * Not a very clean design
     */
    auto const *ev = GetEvent();

    size_t n_pts = 0;
    for (auto const& label : _producer_labels) {
      auto handle = ev->getValidHandle<std::vector<recob::SpacePoint>>(label);

      if (! handle.isValid()) {
          LARCV_ERROR() << "Failed to get SpacePoint from " << label << std::endl;
          return false;
      }

      n_pts += handle->size();
    }

    // reserve
    v_occupancy.reserve(n_pts);
    v_charge.reserve(n_pts);
    v_charge_asym.reserve(n_pts);
    v_chi2.reserve(n_pts);

    if (_store_wire_info) {
      LARCV_DEBUG() << "Store wire info\n";
      for (size_t plane = 0; plane < _n_planes; ++plane) {
        v_hit_charge[plane].reserve(n_pts);
        v_hit_amp[plane].reserve(n_pts);
        v_hit_time[plane].reserve(n_pts);
        v_hit_rms[plane].reserve(n_pts);
      }
    }

    for (auto const& label : _producer_labels) {
      auto handle = ev->getValidHandle<std::vector<recob::SpacePoint>>(label);
      auto const &points = *handle;

      art::InputTag const producer_tag(label);
      art::FindManyP<recob::Hit> find_hits(handle, *ev, producer_tag);

      size_t n_dropped = 0;
      for (size_t i_pt = 0; i_pt < points.size(); ++i_pt) {
          if (i_pt < 10
              || (i_pt < 1000 && i_pt % 100 == 0)
              || (i_pt < 10000 && i_pt % 1000 ==0)
              || (i_pt % 10000 == 0)) {
            LARCV_INFO() << "Processing " << label << ": " 
              << i_pt << "/" << points.size() << '\n';
          }

          auto const &pt = points[i_pt];

          // calculation from Tracys' Cluster3D
          float charge = pt.ErrXYZ()[1];
          float charge_asym = pt.ErrXYZ()[3];

          if (charge < _reco_charge_range[0] || charge > _reco_charge_range[1]) {
            continue;
          }

          auto *xyz = pt.XYZ();
          VoxelID_t vox_id = meta.id(xyz[0], xyz[1], xyz[2]);
          if(vox_id == larcv::kINVALID_VOXELID) { 
              if (n_dropped < _max_debug_dropping)
                  LARCV_DEBUG() << "Dropping space point ("
                      << xyz[0] << ","
                      << xyz[1] << ","
                      << xyz[2] << ")"
                      << std::endl;
              ++n_dropped;
              continue;
          } 

          v_occupancy.emplace(vox_id, 1, true);

          if (!(v_chi2.find(vox_id) == larcv::kINVALID_VOXEL))
              continue;

          v_chi2.emplace(vox_id, pt.Chisq(), true);
          v_charge.emplace(vox_id, charge, true);
          v_charge_asym.emplace(vox_id, charge_asym, true);

          if (_store_wire_info) {
            std::vector<art::Ptr<recob::Hit>> hits;
            if (hits.size() > _n_planes) {
                LARCV_WARNING() 
                    << "Dropping space point - "
                    << "Wrong number of hits: "
                    << hits.size()
                    << " (expecting " << _n_planes << ")" 
                    << std::endl;
                ++n_dropped;
                continue;
            }

            find_hits.get(i_pt, hits);
            for (const auto& hit : hits) {
                size_t plane = hit->WireID().Plane;
                if (plane < 0 || plane >= _n_planes) {
                    LARCV_CRITICAL() << "Invalid plane " << plane << std::endl;
                    continue;
                }

                float charge  = hit->Integral();
                v_hit_charge[plane].emplace(vox_id, charge,               true);
                v_hit_amp   [plane].emplace(vox_id, hit->PeakAmplitude(), true);
                v_hit_time  [plane].emplace(vox_id, hit->PeakTime(),      true);
                v_hit_rms   [plane].emplace(vox_id, hit->RMS(),           true);
            }
         }
      }

      LARCV_INFO() << n_dropped << " out of " << points.size() 
          << " SpacePoints dropped from " << label 
          << std::endl;
    }

    auto store = [&](auto &vset, const std::string& name)
    {
        if (_drop_output.count(name) == 1) return;
        //auto &tensor = reinterpret_cast<larcv::EventSparseTensor3D>(
        std::string label = _output_label;
	      if(!name.empty()) label = label + "_" + name;
        auto &tensor = mgr.get_data<larcv::EventSparseTensor3D>(label);
        tensor.emplace(std::move(vset), meta);
    };

    auto store_vec = [&](auto &vec, const std::string& name)
    {
        if (_drop_output.count(name) == 1) return;
        for (size_t i = 0; i < _n_planes; ++i)
            store(vec[i], name + std::to_string(i));
    };

    store(v_charge,      "");
    store(v_charge_asym, "charge_asym");
    store(v_chi2,        "chi2");
    store(v_occupancy,   "occupancy");

    if (_store_wire_info) {
      store_vec(v_hit_charge, "hit_charge");
      store_vec(v_hit_amp,    "hit_amp");
      store_vec(v_hit_time,   "hit_time");
      store_vec(v_hit_rms,    "hit_rms");
    }

    return true;
  }
      
  void SuperaSpacePoint::finalize()
  {}

  float SuperaSpacePoint::get_common_charge(
      const std::vector<art::Ptr<recob::Hit>>& hits)
  {
      float charge_common = 0.;
      float charge_total  = 0.;
      int idx0(std::numeric_limits<int>::min());
      int idx1(std::numeric_limits<int>::max());

      for (const auto& hit : hits)  {
          charge_total += hit->Integral();

          float peak = hit->PeakTime();
          float width = 2. * hit->RMS() + 0.5;
          int start = peak - width;
          int stop  = peak + width;
          idx0 = std::max(start,    idx0);
          idx1 = std::min(stop + 1, idx1);
      }
      
      /* std::cout << "Idx " << idx0 << " -> " << idx1 << std::endl; */

      auto integrate = [&](const auto& hit)
      {
          double mean  = hit->PeakTime();
          double amp   = hit->PeakAmplitude();
          double width = hit->RMS();

          float t0 = (idx0 - mean) / width / TMath::Sqrt(2.);
          float t1 = (idx1 - mean) / width / TMath::Sqrt(2.);
          float integral = TMath::Erf(t1) - TMath::Erf(t0);
          integral *= amp * TMath::Sqrt(TMath::PiOver2());

           /* for (int pos = idx0; pos < idx1; ++pos)
           integral += amp * TMath::Gaus(double(pos), mean, width); */

          /*
          std::cout
            << "m:" << mean
            << " a:" << amp
            << " w:" << width
            << " q:" << integral
            << std::endl;
            */

          return integral;
      };

      if (charge_total <= 0 || idx0 >= idx1)
        return 0;

      for (const auto& hit : hits)
        charge_common += integrate(hit);

      return charge_common;
  }
}

#endif
