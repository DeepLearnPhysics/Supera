#ifndef __SUPERA_TYPES_H__
#define __SUPERA_TYPES_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
namespace supera {

  /// enum to define LAr* data type
  enum class LArDataType_t : unsigned int {
    kLArWire_t,                 ///< recob::Wire
    kLArHit_t,                  ///< recob::Hit
    kLArOpDigit_t,              ///< raw::OpDetWaveform
    kLArMCTruth_t,              ///< simb::MCTruth
    kLArMCParticle_t,           ///< simb::MCParticle
    kLArMCMiniPart_t,           ///< sim::MCParticleLite
    kLArMCTrack_t,              ///< sim::MCTrack
    kLArMCShower_t,             ///< sim::MCShower
    kLArSimCh_t,                ///< sim::SimChannel
    kLArSimEnergyDeposit_t,     ///< sim::SimEnergyDeposit
    kLArSimEnergyDepositLite_t, ///< sim::SimEnergyDepositLite
    kLArSpacePoint_t,           ///< simb::SpacePoint
    kLArOpFlash_t,              ///< recob::OpFlash
    kLArCRTHit_t,               ///< sbn::crt::CRTHit
    kLArDataTypeMax
  };

  template <class T>
  LArDataType_t LArDataType();

  template <>
  LArDataType_t LArDataType<supera::LArHit_t>();
  template <>
  LArDataType_t LArDataType<supera::LArWire_t>();
  template <>
  LArDataType_t LArDataType<supera::LArOpDigit_t>();
  template <>
  LArDataType_t LArDataType<supera::LArMCTruth_t>();
  template <>
  LArDataType_t LArDataType<supera::LArMCParticle_t>();
  template <>
  LArDataType_t LArDataType<supera::LArMCMiniPart_t>();
  template <>
  LArDataType_t LArDataType<supera::LArMCTrack_t>();
  template <>
  LArDataType_t LArDataType<supera::LArMCShower_t>();
  template <>
  LArDataType_t LArDataType<supera::LArSimCh_t>();
  template <>
  LArDataType_t LArDataType<supera::LArSimEnergyDeposit_t>();
  template <>
  LArDataType_t LArDataType<supera::LArSimEnergyDepositLite_t>();
  template <>
  LArDataType_t LArDataType<supera::LArSpacePoint_t>();
  template <>
  LArDataType_t LArDataType<supera::LArOpFlash_t>();
  template <>
  LArDataType_t LArDataType<supera::LArCRTHit_t>();

  class RSEID {
  public:
    RSEID(size_t run_val = 0, size_t subrun_val = 0, size_t event_val = 0)
      : run(run_val), subrun(subrun_val), event(event_val)
    {}
    ~RSEID() {}

    inline bool operator<(const RSEID& rhs) const
    {
      if (run < rhs.run) return true;
      if (run > rhs.run) return false;
      if (subrun < rhs.subrun) return true;
      if (subrun > rhs.subrun) return false;
      if (event < rhs.event) return true;
      if (event > rhs.event) return false;
      return false;
    }

    inline bool operator==(const RSEID& rhs) const
    {
      return (run == rhs.run && subrun == rhs.subrun && event == rhs.event);
    }

    inline bool operator!=(const RSEID& rhs) const { return !((*this) == rhs); }

    inline bool operator>(const RSEID& rhs) const { return ((*this) != rhs && !((*this) < rhs)); }

    size_t run, subrun, event;
  };

}
//#endif
//#endif
#endif
