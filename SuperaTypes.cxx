#ifndef __SUPERA_TYPES_CXX__
#define __SUPERA_TYPES_CXX__

#include "SuperaTypes.h"
namespace supera {
  template <>
  LArDataType_t LArDataType<supera::LArHit_t>()
  {
    return LArDataType_t::kLArHit_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArWire_t>()
  {
    return LArDataType_t::kLArWire_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArOpDigit_t>()
  {
    return LArDataType_t::kLArOpDigit_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArMCParticle_t>()
  {
    return LArDataType_t::kLArMCParticle_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArMCMiniPart_t>()
  {
    return LArDataType_t::kLArMCMiniPart_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArMCTruth_t>()
  {
    return LArDataType_t::kLArMCTruth_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArMCTrack_t>()
  {
    return LArDataType_t::kLArMCTrack_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArMCShower_t>()
  {
    return LArDataType_t::kLArMCShower_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArSimCh_t>()
  {
    return LArDataType_t::kLArSimCh_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArSimEnergyDeposit_t>()
  {
    return LArDataType_t::kLArSimEnergyDeposit_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArSimEnergyDepositLite_t>()
  {
    return LArDataType_t::kLArSimEnergyDepositLite_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArSpacePoint_t>()
  {
    return LArDataType_t::kLArSpacePoint_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArOpFlash_t>()
  {
    return LArDataType_t::kLArOpFlash_t;
  }
  template <>
  LArDataType_t LArDataType<supera::LArCRTHit_t>()
  {
    return LArDataType_t::kLArCRTHit_t;
  }
}
#endif
