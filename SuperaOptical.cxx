#ifndef __SUPERAOPTICAL_CXX__
#define __SUPERAOPTICAL_CXX__

#include "SuperaOptical.h"
#include "larcv/core/DataFormat/EventFlash.h"

namespace larcv {

  static SuperaOpticalProcessFactory __global_SuperaOpticalProcessFactory__;

  SuperaOptical::SuperaOptical(const std::string name) : SuperaBase(name) {}

  void SuperaOptical::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _opflash_producer_label_v = cfg.get<std::vector<std::string>>("OpFlashProducers");
    _opflash_output_label_v = cfg.get<std::vector<std::string>>("OpFlashOutputs");

    if (_opflash_producer_label_v.size() != _opflash_output_label_v.size()) {
      LARCV_CRITICAL() << "Producer and output labels need to be the same size" << std::endl;
      throw larbys();
    }

    for (auto const& label : _opflash_producer_label_v)
      if (!label.empty()) Request(supera::LArDataType_t::kLArOpFlash_t, label);
  }

  void SuperaOptical::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaOptical::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto const* ev = GetEvent();

    for (size_t label_idx = 0; label_idx < _opflash_output_label_v.size(); ++label_idx) {
      std::string _opflash_output_label = _opflash_output_label_v[label_idx];
      std::string _opflash_producer_label = _opflash_producer_label_v[label_idx];

      auto& opflash_tensor = mgr.get_data<larcv::EventFlash>(_opflash_output_label);
      //auto const& meta = opflash_tensor.meta();

      auto handle = ev->getValidHandle<std::vector<supera::LArOpFlash_t>>(_opflash_producer_label);
      if (!handle.isValid()) {
        LARCV_WARNING() << "No valid OpFlash found for label " << _opflash_producer_label
                        << std::endl;
        return false;
      }

      auto const& flashes = *handle;

      std::vector<larcv::Flash> fset;
      for (size_t idx = 0; idx < flashes.size(); ++idx) {
        auto const& flash = flashes[idx];

        larcv::Flash larcv_flash(flash.Time(),
                                 flash.TimeWidth(),
                                 flash.AbsTime(),
                                 flash.Frame(),
                                 flash.PEs(),
                                 flash.InBeamFrame(),
                                 flash.OnBeamTime(),
                                 flash.FastToTotal(),
                                 flash.XCenter(),
                                 flash.XWidth(),
                                 flash.YCenter(),
                                 flash.YWidth(),
                                 flash.ZCenter(),
                                 flash.ZWidth(),
                                 flash.WireCenters(),
                                 flash.WireWidths(),
                                 idx);
        fset.push_back(larcv_flash);
      }

      opflash_tensor.emplace(std::move(fset));
    }
    return true;
  }

  void SuperaOptical::finalize() {}
}
#endif
