#ifndef __SUPERABBOXINTERACTION_CXX__
#define __SUPERABBOXINTERACTION_CXX__

#include "SuperaBBoxInteraction.h"
#include "GenRandom.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaBBoxInteractionProcessFactory __global_SuperaBBoxInteractionProcessFactory__;

  SuperaBBoxInteraction::SuperaBBoxInteraction(const std::string name) : SuperaBase(name) {}
  void SuperaBBoxInteraction::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _use_sed_lite = cfg.get<bool>("UseSEDLite", false);

    _origin = cfg.get<unsigned short>("Origin", 0);

    _cluster3d_labels = cfg.get<std::vector<std::string>>("Cluster3DLabels");
    _tensor3d_labels = cfg.get<std::vector<std::string>>("Tensor3DLabels");

    auto bbox_size = cfg.get<std::vector<double>>("BBoxSize");
    assert(bbox_size.size() == 3);
    _xlen = bbox_size.at(0);
    _ylen = bbox_size.at(1);
    _zlen = bbox_size.at(2);

    auto voxel_size = cfg.get<std::vector<double>>("VoxelSize");
    assert(voxel_size.size() == 3);
    _xvox = voxel_size.at(0);
    _yvox = voxel_size.at(1);
    _zvox = voxel_size.at(2);

    _use_fixed_bbox = cfg.get<bool>("UseFixedBBox", false);
    _bbox_bottom = cfg.get<std::vector<double>>("BBoxBottom", {0, 0, 0});
    assert(_bbox_bottom.size() == 3);

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

  void SuperaBBoxInteraction::initialize()
  {
    SuperaBase::initialize();
  }

  template <class T>
  void SuperaBBoxInteraction::AdaptBBoxToSED(std::vector<T> const& sedep_v, larcv::BBox3D& bbox)
  {
    LARCV_INFO() << "Processing SimEnergyDeposit array: " << sedep_v.size() << std::endl;
    for (size_t sedep_idx = 0; sedep_idx < sedep_v.size(); ++sedep_idx) {
      auto const& sedep = sedep_v.at(sedep_idx);
      larcv::Point3D pt;
      pt.x = sedep.X();
      pt.y = sedep.Y();
      pt.z = sedep.Z();
      if (!update_bbox(bbox, pt)) break;
    }
  }

  bool SuperaBBoxInteraction::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);
    larcv::BBox3D bbox(0, 0, 0, 0, 0, 0);

    if (_use_fixed_bbox) {
      double x0 = _bbox_bottom[0];
      double y0 = _bbox_bottom[1];
      double z0 = _bbox_bottom[2];

      double x1 = x0 + _xlen;
      double y1 = y0 + _ylen;
      double z1 = z0 + _zlen;

      larcv::Point3D p0(x0, y0, z0);
      larcv::Point3D p1(x1, y1, z1);
      bbox.update(p0, p1);
    }
    else {
      /*
    // Retrieve mcparticles
    art::Handle<std::vector<simb::MCParticle> > mcpHandle;
    evt.getByLabel(fMCParticleLabel,mcpHandle);
    if(!mcpHandle.isValid()) throw cet::exception(__FUNCTION__) << "Failed to retrieve simb::MCParticle";;

    // Find associations
    art::FindOneP<simb::MCTruth> ass(mcpHandle, evt, fMCParticleLabel);
    std::vector<simb::Origin_t> orig_array;
    orig_array.reserve(mcpHandle->size());
    for(size_t i=0; i<mcpHandle->size(); ++i) {
      const art::Ptr<simb::MCTruth> &mct = ass.at(i);
      orig_array.push_back(mct->Origin());
    }

    const std::vector<simb::MCParticle>& mcp_array(*mcpHandle);
    fPart.AddParticles(mcp_array,orig_array);
    */
      // Register primary vertex points
      auto const& mct_v = LArData<supera::LArMCTruth_t>();
      LARCV_INFO() << "Processing MCTruth: Loaded MCTruth array: " << mct_v.size() << std::endl;
      for (size_t mct_index = 0; mct_index < mct_v.size(); ++mct_index) {
        auto const& mct = mct_v[mct_index];
        if (_origin && mct.Origin() != _origin) {
          LARCV_INFO() << "Skipping MCTruth of oritin type: " << (int)(mct.Origin()) << std::endl;
          continue;
        }

        for (int i = 0; i < mct.NParticles(); ++i) {
          auto const& mcp = mct.GetParticle(i);
          if (mcp.StatusCode() != 1) {
            LARCV_INFO() << "Skipping MCTruth::MCParticle of status code: "
                         << (int)(mcp.StatusCode()) << std::endl;
            continue;
          }

          auto const& pos = mcp.Position(0);
          larcv::Point3D pt;
          pt.x = pos.X();
          pt.y = pos.Y();
          pt.z = pos.Z();
          LARCV_INFO() << "Registering MCTruth::MCParticle vertex: (" << pt.x << "," << pt.y << ","
                       << pt.z << ")"
                       << " ... PDG " << mcp.PdgCode() << std::endl;
          this->update_bbox(bbox, pt);
        }
      }

      // Register particle energy deposition coordinates
      if (_use_sed_lite) {
        AdaptBBoxToSED<sim::SimEnergyDepositLite>(LArData<supera::LArSimEnergyDepositLite_t>(),
                                                  bbox);
      }
      else {
        AdaptBBoxToSED<sim::SimEnergyDeposit>(LArData<supera::LArSimEnergyDeposit_t>(), bbox);
      }

      // Randomize BBox location
      randomize_bbox_center(bbox);
    }

    // Create 3D meta
    larcv::Voxel3DMeta meta;
    auto const& min_pt = bbox.bottom_left();
    auto const& max_pt = bbox.top_right();
    size_t xnum = _xlen / _xvox;
    size_t ynum = _ylen / _yvox;
    size_t znum = _zlen / _zvox;
    meta.set(min_pt.x, min_pt.y, min_pt.z, max_pt.x, max_pt.y, max_pt.z, xnum, ynum, znum);
    LARCV_INFO() << "3D Meta:" << meta.dump() << std::endl;

    // Create Cluster3D
    for (auto const& name : _cluster3d_labels) {
      auto& cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(name);
      cluster3d.meta(meta);
    }
    // Create Tensor3D
    for (auto const& name : _tensor3d_labels) {
      auto& tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(name);
      tensor3d.meta(meta);
    }
    return true;
  }

  void SuperaBBoxInteraction::finalize() {}

  bool SuperaBBoxInteraction::update_bbox(larcv::BBox3D& bbox, const larcv::Point3D& pt)
  {
    larcv::Point3D min_pt, max_pt;
    if (bbox.empty()) {
      min_pt = max_pt = pt;
      max_pt.x += 1.e-9;
      max_pt.y += 1.e-9;
      max_pt.z += 1.e-9;
      bbox.update(min_pt, max_pt);
      LARCV_INFO() << "Defining minimal BBox:" << bbox.dump();
      return true;
    }
    else {
      min_pt = bbox.bottom_left();
      max_pt = bbox.top_right();
    }
    if (_world_bounds.contains(pt)) {
      if (!bbox.contains(pt)) {
        LARCV_DEBUG() << "Updating BBox:" << bbox.dump();
        if ((max_pt.x - min_pt.x) < _xlen) {
          if (pt.x < min_pt.x) {
            if ((max_pt.x - pt.x) < _xlen)
              min_pt.x = pt.x;
            else
              min_pt.x = max_pt.x - _xlen;
          }
          else if (pt.x > max_pt.x) {
            if ((pt.x - min_pt.x) < _xlen)
              max_pt.x = pt.x;
            else
              max_pt.x = min_pt.x + _xlen;
          }
        }
        if ((max_pt.y - min_pt.y) < _ylen) {
          if (pt.y < min_pt.y) {
            if ((max_pt.y - pt.y) < _ylen)
              min_pt.y = pt.y;
            else
              min_pt.y = max_pt.y - _ylen;
          }
          else if (pt.y > max_pt.y) {
            if ((pt.y - min_pt.y) < _ylen)
              max_pt.y = pt.y;
            else
              max_pt.y = min_pt.y + _ylen;
          }
        }
        if ((max_pt.z - min_pt.z) < _zlen) {
          if (pt.z < min_pt.z) {
            if ((max_pt.z - pt.z) < _zlen)
              min_pt.z = pt.z;
            else
              min_pt.z = max_pt.z - _zlen;
          }
          else if (pt.z > max_pt.z) {
            if ((pt.z - min_pt.z) < _zlen)
              max_pt.z = pt.z;
            else
              max_pt.z = min_pt.z + _zlen;
          }
        }
        // update bbox
        bbox.update(min_pt, max_pt);
        LARCV_DEBUG() << " ... to:" << bbox.dump();
      }
      else {
        LARCV_DEBUG() << "No update in BBox: point already contained!" << std::endl;
      }
    }
    else {
      LARCV_DEBUG() << "No update in BBox: point outside the world boundary" << std::endl;
    }
    return ((max_pt.x - min_pt.x) < _xlen || (max_pt.y - min_pt.y) < _ylen ||
            (max_pt.z - min_pt.z) < _zlen);
  }

  void SuperaBBoxInteraction::randomize_bbox_center(larcv::BBox3D& bbox)
  {

    larcv::Point3D min_pt = bbox.bottom_left();
    larcv::Point3D max_pt = bbox.top_right();
    LARCV_INFO() << "Randomize before:" << bbox.dump() << std::endl;
    // see if box location can be randomized
    if ((max_pt.x - min_pt.x) < _xlen) {
      double xshift = supera::GenRandom::get().Flat(0, _xlen - (max_pt.x - min_pt.x));
      xshift *= (supera::GenRandom::get().Flat(-1., 1.) > 0. ? 1. : -1.);
      if (xshift > 0) {
        max_pt.x += xshift;
        min_pt.x = max_pt.x - _xlen;
      }
      else {
        min_pt.x += xshift;
        max_pt.x = min_pt.x + _xlen;
      }
    }
    if ((max_pt.y - min_pt.y) < _ylen) {
      double yshift = supera::GenRandom::get().Flat(0, _ylen - (max_pt.y - min_pt.y));
      yshift *= (supera::GenRandom::get().Flat(-1., 1.) > 0. ? 1. : -1.);
      if (yshift > 0) {
        max_pt.y += yshift;
        min_pt.y = max_pt.y - _xlen;
      }
      else {
        min_pt.y += yshift;
        max_pt.y = min_pt.y + _xlen;
      }
    }
    if ((max_pt.z - min_pt.z) < _zlen) {
      double zshift = supera::GenRandom::get().Flat(0, _zlen - (max_pt.z - min_pt.z));
      zshift *= (supera::GenRandom::get().Flat(-1., 1.) > 0. ? 1. : -1.);
      if (zshift > 0) {
        max_pt.z += zshift;
        min_pt.z = max_pt.z - _xlen;
      }
      else {
        min_pt.z += zshift;
        max_pt.z = min_pt.z + _xlen;
      }
    }
    bbox.update(min_pt, max_pt);
    LARCV_INFO() << "Randomize after:" << bbox.dump() << std::endl;
  }

}

#endif
