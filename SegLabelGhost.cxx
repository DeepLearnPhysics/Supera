#ifndef __SEGLABELGHOST_CXX__
#define __SEGLABELGHOST_CXX__

#include "SegLabelGhost.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SegLabelGhostProcessFactory __global_SegLabelGhostProcessFactory__;

  SegLabelGhost::SegLabelGhost(const std::string name) : ProcessBase(name) {}

  void SegLabelGhost::configure(const PSet& cfg)
  {
    _reco_prod = cfg.get<std::string>("RecoProducer");
    _mc_prod = cfg.get<std::string>("McProducer");
    _out_prod = cfg.get<std::string>("OutputProducer");
    _min_num_voxel = cfg.get<size_t>("MinVoxelCount", 0);
  }

  void SegLabelGhost::initialize() {}

  bool SegLabelGhost::process(IOManager& mgr)
  {
    auto const& evt_reco = mgr.get_data<larcv::EventSparseTensor3D>(_reco_prod);
    auto const& evt_mc = mgr.get_data<larcv::EventSparseTensor3D>(_mc_prod);

    auto evt_out =
      reinterpret_cast<larcv::EventSparseTensor3D*>(mgr.get_data("sparse3d", _out_prod));

    evt_out->meta(evt_reco.meta());

    for (auto const& reco : evt_reco.as_vector()) {
      auto const& mc = evt_mc.find(reco.id());
      float class_def = mc.id() == larcv::kINVALID_VOXELID ? 1 : 0;
      ((VoxelSet*)evt_out)->emplace(reco.id(), class_def, false);
    }

    if (_min_num_voxel < 1) return true;

    if (evt_out->as_vector().size() < _min_num_voxel) {
      LARCV_NORMAL() << "Skipping event " << evt_out->event_key() << " due to voxel count ("
                     << evt_out->as_vector().size() << " < " << _min_num_voxel << ")" << std::endl;
      return false;
    }
    return true;
  }

  void SegLabelGhost::finalize() {}

}
#endif
