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
    _producer_label     = cfg.get<std::string>("SpacePointProducer");
    _output_label       = cfg.get<std::string>("OutputLabel");
    _max_debug_dropping = cfg.get<size_t>("MaxDebugForDropping", 0);

    Request(supera::LArDataType_t::kLArSpacePoint_t, _producer_label);
  }

  void SuperaSpacePoint::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaSpacePoint::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    // Main cluster3d to be filled
    auto& event_cluster_v = mgr.get_data<larcv::EventClusterVoxel3D>(_output_label);
    auto const& meta = event_cluster_v.meta();
    LARCV_INFO() << "Voxel3DMeta: " << meta.dump();

    std::vector<larcv::VoxelSet> v_charge(1);
    std::vector<larcv::VoxelSet> v_chi2(1);


    /* FIXME(kvtsang) To be removed?
     * Find associated hits
     * Not a very clean design
     */
    auto const *ev = GetEvent();
    auto space_point_handle = 
        ev->getValidHandle<std::vector<recob::SpacePoint>>(_producer_label);

    if (! space_point_handle.isValid()) {
        LARCV_ERROR() << "Failed to get SpacePoint from " 
            << _producer_label << std::endl;
        return false;
    }
    auto const &points = *space_point_handle;
    art::InputTag const producer_tag(_producer_label);
    art::FindManyP<recob::Hit> find_hits(space_point_handle, *ev, producer_tag);

    //TODO(kvtsang) Preferred methond to get data vector via Supera FW
    //auto const& points = LArData<supera::LArSpacePoint_t>();


    size_t n_dropped = 0;
    for (size_t i_pt = 0; i_pt < points.size(); ++i_pt) {
        auto const &pt = points[i_pt];
    //for (auto const &pt : points) {
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
        } 
        else {
            v_chi2[0].emplace(vox_id, pt.Chisq(), true);

            /* Calculuate charge by arverage 3 wires
             * FIXME(kvtsang) should be provided by SpacePoint
             */ 
            std::vector<art::Ptr<recob::Hit>> hits;
            find_hits.get(i_pt, hits);
            float charge = average_hit_charge(hits);
            v_charge[0].emplace(vox_id, charge, true);
        }
    }

    LARCV_INFO() << n_dropped << " out of " << points.size() 
        << " SpacePoints dropped."
        << std::endl;

    auto store = [&](auto &vec, const std::string& name) {
        auto &cluster = 
            mgr.get_data<larcv::EventClusterVoxel3D>(_output_label + name);
        larcv::VoxelSetArray vsa;
        vsa.emplace(std::move(vec));
        cluster.emplace(std::move(vsa), meta);
    };

    store(v_charge, ""     );
    store(v_chi2,   "_chi2");

    return true;
  }
      
  void SuperaSpacePoint::finalize()
  {}

  float SuperaSpacePoint::average_hit_charge(
          const std::vector<art::Ptr<recob::Hit>>& hits)
  {
      float charge(0.);
      int idx0(std::numeric_limits<int>::min());
      int idx1(std::numeric_limits<int>::max());

      for (const auto& hit : hits)  {
          charge += hit->Integral();

          float peak = hit->PeakTime();
          float width = 2. * hit->RMS() + 0.5;
          int start = peak - width;
          int stop  = peak + width;
          idx0 = std::max(start,    idx0);
          idx1 = std::min(stop + 1, idx1);
      }

      if (!hits.empty())
          charge /= float(hits.size());

      auto integrate = [&](const auto& hit) {
          double mean  = hit->PeakTime();
          double amp   = hit->PeakAmplitude();
          double width = hit->RMS();

          double integral(0.);
          for (int pos = idx0; pos < idx1; ++pos)
              integral += amp * TMath::Gaus(double(pos), mean, width);
          return integral;
      };

      float integral(0.);
      if (charge > 0 && idx1 > idx0) {
          for (const auto& hit : hits) {
              integral += integrate(hit);
          }
      }
      return integral;
  }
}

#endif
