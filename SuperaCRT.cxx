#ifndef __SUPERACRT_CXX__
#define __SUPERACRT_CXX__

#include "SuperaCRT.h"
#include "larcv/core/DataFormat/EventCRTHit.h"

namespace larcv {

  static SuperaCRTProcessFactory __global_SuperaCRTProcessFactory__;

  SuperaCRT::SuperaCRT(const std::string name) : SuperaBase(name) {}

  void SuperaCRT::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _crthit_producer_label_v = cfg.get<std::vector<std::string>>("CRTHitProducers");
    _crthit_output_label_v = cfg.get<std::vector<std::string>>("CRTHitOutputs");
    _use_crtsp = cfg.get<bool>("UseCRTSpacePoint", false);

    if (_crthit_producer_label_v.size() != _crthit_output_label_v.size()) {
      LARCV_CRITICAL() << "Producer and output labels need to be the same size" << std::endl;
      throw larbys();
    }

    for (auto const& label : _crthit_producer_label_v)
      if (!label.empty()) {
        if (!_use_crtsp) {
          Request(supera::LArDataType_t::kLArCRTHit_t, label);
        } else {
          Request(supera::LArDataType_t::kLArCRTSpacePoint_t, label);
        }
      }
  }

  void SuperaCRT::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaCRT::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto const* ev = GetEvent();

    for (size_t label_idx = 0; label_idx < _crthit_output_label_v.size(); ++label_idx) {
      std::string _crthit_output_label = _crthit_output_label_v[label_idx];
      std::string _crthit_producer_label = _crthit_producer_label_v[label_idx];

      auto& crthit_tensor = mgr.get_data<larcv::EventCRTHit>(_crthit_output_label);
      //auto const& meta = opflash_tensor.meta();

      // Declare a variable to hold the data
      bool valid_handle = false;
      std::vector<larcv::CRTHit> cset;
      
      if (_use_crtsp) {
        auto crt_sp_handle = ev->getValidHandle<std::vector<supera::LArCRTSpacePoint_t>>(_crthit_producer_label);
        valid_handle = true;  // getValidHandle throws if not valid, so we'll never get here if invalid
        
        auto const& crthits = *crt_sp_handle;
        for (size_t idx = 0; idx < crthits.size(); ++idx) {
          auto const& crthit = crthits[idx];
          larcv::CRTHit larcv_crthit;
          larcv_crthit.id(idx);
          
          larcv_crthit.peshit(crthit.PE());
          larcv_crthit.ts0_s(0);
          larcv_crthit.ts0_ns(crthit.Ts0());
          larcv_crthit.ts1_ns(crthit.Ts1());
          larcv_crthit.x_pos(crthit.X());
          larcv_crthit.y_pos(crthit.Y());
          larcv_crthit.z_pos(crthit.Z());
          larcv_crthit.x_err(crthit.XErr());
          larcv_crthit.y_err(crthit.YErr());
          larcv_crthit.z_err(crthit.ZErr());
          
          cset.push_back(larcv_crthit);
        }
      } else {
        auto crt_hit_handle = ev->getValidHandle<std::vector<supera::LArCRTHit_t>>(_crthit_producer_label);
        valid_handle = true;
        
        auto const& crthits = *crt_hit_handle;
        for (size_t idx = 0; idx < crthits.size(); ++idx) {
          auto const& crthit = crthits[idx];
          larcv::CRTHit larcv_crthit;
          larcv_crthit.id(idx);
          
          larcv_crthit.feb_id(crthit.feb_id);
          larcv_crthit.pesmap(crthit.pesmap);
          larcv_crthit.peshit(crthit.peshit);
          larcv_crthit.ts0_s(crthit.ts0_s);
          larcv_crthit.ts0_s_corr(crthit.ts0_s_corr);
          larcv_crthit.ts0_ns(crthit.ts0_ns);
          larcv_crthit.ts0_ns_corr(crthit.ts0_ns_corr);
          larcv_crthit.ts1_ns(crthit.ts1_ns);
          larcv_crthit.plane(crthit.plane);
          larcv_crthit.x_pos(crthit.x_pos);
          larcv_crthit.x_err(crthit.x_err);
          larcv_crthit.y_pos(crthit.y_pos);
          larcv_crthit.y_err(crthit.y_err);
          larcv_crthit.z_pos(crthit.z_pos);
          larcv_crthit.z_err(crthit.z_err);
          larcv_crthit.tagger(crthit.tagger);
          
          cset.push_back(larcv_crthit);
        }
      }

      // Handle the case where neither was valid (though this shouldn't happen with getValidHandle)
      if (!valid_handle) {
        LARCV_WARNING() << "No valid CRT data found for label " << _crthit_producer_label << std::endl;
        return false;
      }

      crthit_tensor.emplace(std::move(cset));
    }
    return true;
  }

  void SuperaCRT::finalize() {}
}

#endif
