#ifndef __SUPERABASE_CXX__
#define __SUPERABASE_CXX__

#include "SuperaBase.h"

namespace larcv {

  static SuperaBaseProcessFactory __global_SuperaBaseProcessFactory__;

  SuperaBase::SuperaBase(const std::string name) : ProcessBase(name), _empty_string()
  {
    ClearEventData();
  }

  void SuperaBase::Request(supera::LArDataType_t type, std::string name)
  {
    _data_request_m[type] = name;
  }

  const std::string& SuperaBase::LArDataLabel(supera::LArDataType_t type) const
  {
    auto iter = _data_request_m.find(type);
    if (iter == _data_request_m.end()) return _empty_string;
    return (*iter).second;
  }

  void SuperaBase::configure(const PSet& cfg)
  {
    _time_offset = cfg.get<int>("TimeOffset", 2400);

    auto producer_wire = cfg.get<std::string>("LArWireProducer", "");
    auto producer_hit = cfg.get<std::string>("LArHitProducer", "");
    auto producer_opdigit = cfg.get<std::string>("LArOpDigitProducer", "");
    auto producer_mctruth = cfg.get<std::string>("LArMCTruthProducer", "");
    auto producer_mcpart = cfg.get<std::string>("LArMCParticleProducer", "");
    auto producer_mcmp = cfg.get<std::string>("LArMCMiniPartProducer", "");
    auto producer_mctrack = cfg.get<std::string>("LArMCTrackProducer", "");
    auto producer_mcshower = cfg.get<std::string>("LArMCShowerProducer", "");
    auto producer_simch = cfg.get<std::string>("LArSimChProducer", "");
    auto producer_simedep = cfg.get<std::string>("LArSimEnergyDepositProducer", "");
    auto producer_simedep_lite = cfg.get<std::string>("LArSimEnergyDepositLiteProducer", "");
    auto producer_sps = cfg.get<std::string>("LArSpacePoint", "");
    auto producer_opflash = cfg.get<std::string>("LArOpFlashProducer", "");
    auto producer_crthit = cfg.get<std::string>("LArCRTHitProducer", "");

    if (!producer_wire.empty()) {
      LARCV_INFO() << "Requesting Wire data product by " << producer_wire << std::endl;
      Request(supera::LArDataType_t::kLArWire_t, producer_wire);
    }

    if (!producer_hit.empty()) {
      LARCV_INFO() << "Requesting Hit data product by " << producer_hit << std::endl;
      Request(supera::LArDataType_t::kLArHit_t, producer_hit);
    }

    if (!producer_opdigit.empty()) {
      LARCV_INFO() << "Requesting OpDigit data product by " << producer_opdigit << std::endl;
      Request(supera::LArDataType_t::kLArOpDigit_t, producer_opdigit);
    }

    if (!producer_mctruth.empty()) {
      LARCV_INFO() << "Requesting MCTruth data product by " << producer_mctruth << std::endl;
      Request(supera::LArDataType_t::kLArMCTruth_t, producer_mctruth);
    }

    if (!producer_mcpart.empty()) {
      LARCV_INFO() << "Requesting MCParticle data product by " << producer_mcpart << std::endl;
      Request(supera::LArDataType_t::kLArMCParticle_t, producer_mcpart);
    }

    if (!producer_mcmp.empty()) {
      LARCV_INFO() << "Requesting MCMiniPart data product by " << producer_mcmp << std::endl;
      Request(supera::LArDataType_t::kLArMCMiniPart_t, producer_mcmp);
    }

    if (!producer_mctrack.empty()) {
      LARCV_INFO() << "Requesting MCTrack data product by " << producer_mctrack << std::endl;
      Request(supera::LArDataType_t::kLArMCTrack_t, producer_mctrack);
    }

    if (!producer_mcshower.empty()) {
      LARCV_INFO() << "Requesting MCShower data product by " << producer_mcshower << std::endl;
      Request(supera::LArDataType_t::kLArMCShower_t, producer_mcshower);
    }

    if (!producer_simch.empty()) {
      LARCV_INFO() << "Requesting SimCh data product by " << producer_simch << std::endl;
      Request(supera::LArDataType_t::kLArSimCh_t, producer_simch);
    }

    if (!producer_simedep.empty()) {
      LARCV_INFO() << "Requesting SimEnergyDeposit data product by " << producer_simedep
                   << std::endl;
      Request(supera::LArDataType_t::kLArSimEnergyDeposit_t, producer_simedep);
    }

    if (!producer_simedep_lite.empty()) {
      LARCV_INFO() << "Requesting SimEnergyDepositLite data product by " << producer_simedep_lite
                   << std::endl;
      Request(supera::LArDataType_t::kLArSimEnergyDepositLite_t, producer_simedep_lite);
    }

    if (!producer_sps.empty()) {
      LARCV_INFO() << "Requesting SpacePoint data product by " << producer_sps << std::endl;
      Request(supera::LArDataType_t::kLArSpacePoint_t, producer_sps);
    }

    if (!producer_opflash.empty()) {
      LARCV_INFO() << "Requesting Opflash data product by " << producer_opflash << std::endl;
      Request(supera::LArDataType_t::kLArOpFlash_t, producer_opflash);
    }

    if (!producer_crthit.empty()) {
      LARCV_INFO() << "Requesting CRTHit data product by " << producer_crthit << std::endl;
      Request(supera::LArDataType_t::kLArCRTHit_t, producer_crthit);
    }
  }

  void SuperaBase::initialize()
  {
    ClearEventData();
  }

  bool SuperaBase::process(IOManager& mgr)
  {
    return true;
  }

  void SuperaBase::finalize()
  {
    ClearEventData();
  }

  bool SuperaBase::is(const std::string question) const
  {
    if (question == "Supera") return true;
    return false;
  }

  void SuperaBase::ClearEventData()
  {
    _ptr_wire_v = nullptr;
    _ptr_hit_v = nullptr;
    _ptr_opdigit_v = nullptr;
    _ptr_sch_v = nullptr;
    _ptr_mctruth_v = nullptr;
    _ptr_mcp_v = nullptr;
    _ptr_mcmp_v = nullptr;
    _ptr_mct_v = nullptr;
    _ptr_mcs_v = nullptr;
    _ptr_simedep_v = nullptr;
    _ptr_simedep_lite_v = nullptr;
    _ptr_opflash_v = nullptr;
    _ptr_crthit_v = nullptr;

    // FIXME(kvtsang) Temporary solution to access associations
    _event = nullptr;
  }

  // FIXME(kvtsang) Temporary solution to access associations
  const art::Event* SuperaBase::GetEvent()
  {
    if (!_event) throw larbys("art::Event not set!");
    return _event;
  }

  template <>
  const std::vector<supera::LArWire_t>& SuperaBase::LArData<supera::LArWire_t>() const
  {
    if (!_ptr_wire_v) throw larbys("Wire data pointer not available");
    return *_ptr_wire_v;
  }

  template <>
  const std::vector<supera::LArHit_t>& SuperaBase::LArData<supera::LArHit_t>() const
  {
    if (!_ptr_hit_v) throw larbys("Hit data pointer not available");
    return *_ptr_hit_v;
  }

  template <>
  const std::vector<supera::LArOpDigit_t>& SuperaBase::LArData<supera::LArOpDigit_t>() const
  {
    if (!_ptr_opdigit_v) throw larbys("OpDigit data pointer not available");
    return *_ptr_opdigit_v;
  }

  template <>
  const std::vector<supera::LArMCTruth_t>& SuperaBase::LArData<supera::LArMCTruth_t>() const
  {
    if (!_ptr_mctruth_v) throw larbys("MCTruth data pointer not available");
    return *_ptr_mctruth_v;
  }

  template <>
  const std::vector<supera::LArMCParticle_t>& SuperaBase::LArData<supera::LArMCParticle_t>() const
  {
    if (!_ptr_mcp_v) throw larbys("MCParticle data pointer not available");
    return *_ptr_mcp_v;
  }

  template <>
  const std::vector<supera::LArMCMiniPart_t>& SuperaBase::LArData<supera::LArMCMiniPart_t>() const
  {
    if (!_ptr_mcmp_v) throw larbys("MCMiniPart data pointer not available");
    return *_ptr_mcmp_v;
  }

  template <>
  const std::vector<supera::LArMCTrack_t>& SuperaBase::LArData<supera::LArMCTrack_t>() const
  {
    if (!_ptr_mct_v) throw larbys("MCTrack data pointer not available");
    return *_ptr_mct_v;
  }

  template <>
  const std::vector<supera::LArMCShower_t>& SuperaBase::LArData<supera::LArMCShower_t>() const
  {
    if (!_ptr_mcs_v) throw larbys("MCShower data pointer not available");
    return *_ptr_mcs_v;
  }

  template <>
  const std::vector<supera::LArSimCh_t>& SuperaBase::LArData<supera::LArSimCh_t>() const
  {
    if (!_ptr_sch_v) throw larbys("SimCh data pointer not available");
    return *_ptr_sch_v;
  }

  template <>
  const std::vector<supera::LArSimEnergyDeposit_t>&
  SuperaBase::LArData<supera::LArSimEnergyDeposit_t>() const
  {
    if (!_ptr_simedep_v) throw larbys("SimEnergyDeposit data pointer not available");
    return *_ptr_simedep_v;
  }

  template <>
  const std::vector<supera::LArSimEnergyDepositLite_t>&
  SuperaBase::LArData<supera::LArSimEnergyDepositLite_t>() const
  {
    if (!_ptr_simedep_lite_v) throw larbys("SimEnergyDepositLite data pointer not available");
    return *_ptr_simedep_lite_v;
  }

  template <>
  const std::vector<supera::LArSpacePoint_t>& SuperaBase::LArData<supera::LArSpacePoint_t>() const
  {
    if (!_ptr_spacepoint_v) throw larbys("SpacePoint data pointer not available");
    return *_ptr_spacepoint_v;
  }

  template <>
  const std::vector<supera::LArOpFlash_t>& SuperaBase::LArData<supera::LArOpFlash_t>() const
  {
    if (!_ptr_opflash_v) throw larbys("OpFlash data pointer not available");
    return *_ptr_opflash_v;
  }

  template <>
  const std::vector<supera::LArCRTHit_t>& SuperaBase::LArData<supera::LArCRTHit_t>() const
  {
    if (!_ptr_crthit_v) throw larbys("CRTHit data pointer not available");
    return *_ptr_crthit_v;
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArWire_t>& data_v)
  {
    _ptr_wire_v = (std::vector<supera::LArWire_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArHit_t>& data_v)
  {
    _ptr_hit_v = (std::vector<supera::LArHit_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArOpDigit_t>& data_v)
  {
    _ptr_opdigit_v = (std::vector<supera::LArOpDigit_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCTruth_t>& data_v)
  {
    _ptr_mctruth_v = (std::vector<supera::LArMCTruth_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCParticle_t>& data_v)
  {
    _ptr_mcp_v = (std::vector<supera::LArMCParticle_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCMiniPart_t>& data_v)
  {
    _ptr_mcmp_v = (std::vector<supera::LArMCMiniPart_t>*)(&data_v);
    // This assumes that we only load LArMCMiniPart once.
    // Note that these pointers are not reset at each event as of now.
    if (_ptr_mcp_v) {
      for (auto const& mcmp : (*_ptr_mcmp_v))
        _ptr_mcp_v->push_back(supera::LArMCParticle_t(mcmp));
    }
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCTrack_t>& data_v)
  {
    _ptr_mct_v = (std::vector<supera::LArMCTrack_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCShower_t>& data_v)
  {
    _ptr_mcs_v = (std::vector<supera::LArMCShower_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSimCh_t>& data_v)
  {
    _ptr_sch_v = (std::vector<supera::LArSimCh_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSimEnergyDeposit_t>& data_v)
  {
    _ptr_simedep_v = (std::vector<supera::LArSimEnergyDeposit_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSimEnergyDepositLite_t>& data_v)
  {
    _ptr_simedep_lite_v = (std::vector<supera::LArSimEnergyDepositLite_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSpacePoint_t>& data_v)
  {
    _ptr_spacepoint_v = (std::vector<supera::LArSpacePoint_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArOpFlash_t>& data_v)
  {
    _ptr_opflash_v = (std::vector<supera::LArOpFlash_t>*)(&data_v);
  }

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArCRTHit_t>& data_v)
  {
    _ptr_crthit_v = (std::vector<supera::LArCRTHit_t>*)(&data_v);
  }
}

#endif
