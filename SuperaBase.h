/**
 * \file SuperaBase.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaBase
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERABASE_H__
#define __SUPERABASE_H__

#include "FMWKInterface.h"
#include "SuperaTypes.h"
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
//#include "ImageMetaMakerBase.h"

// FIXME(kvtsang) Temporary solution to access associations
#include "art/Framework/Principal/Event.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaBase ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaBase : public ProcessBase {

  public:
    /// Default constructor
    SuperaBase(const std::string name = "SuperaBase");

    /// Default destructor
    ~SuperaBase() {}

    virtual void configure(const PSet&);

    virtual void initialize();

    virtual bool process(IOManager& mgr);

    virtual void finalize();

    virtual bool is(const std::string question) const;

    void Request(supera::LArDataType_t, std::string);

    template <class T>
    void LArData(const T&);

    inline void SetCSV(const std::string& fname) { _csv_fname = fname; }

    // FIXME(kvtsang) Temporary solution to access associations
    void SetEvent(const art::Event* ev) { _event = ev; };

    void ClearEventData();

    //
    // Getter
    //
    const std::string& LArDataLabel(supera::LArDataType_t type) const;

    template <class T>
    const std::vector<T>& LArData() const;

    int TimeOffset() const { return _time_offset; }

    const std::string& CSV() const { return _csv_fname; }

    const art::Event* GetEvent();

  private:
    std::string _empty_string;
    int _time_offset;
    std::map<supera::LArDataType_t, std::string> _data_request_m;

    std::vector<supera::LArWire_t>* _ptr_wire_v;
    std::vector<supera::LArHit_t>* _ptr_hit_v;
    std::vector<supera::LArOpDigit_t>* _ptr_opdigit_v;
    std::vector<supera::LArSimCh_t>* _ptr_sch_v;
    std::vector<supera::LArMCTruth_t>* _ptr_mctruth_v;
    std::vector<supera::LArMCParticle_t>* _ptr_mcp_v;
    std::vector<supera::LArMCMiniPart_t>* _ptr_mcmp_v;
    std::vector<supera::LArMCTrack_t>* _ptr_mct_v;
    std::vector<supera::LArMCShower_t>* _ptr_mcs_v;
    std::vector<supera::LArSimEnergyDeposit_t>* _ptr_simedep_v;
    std::vector<supera::LArSimEnergyDepositLite_t>* _ptr_simedep_lite_v;
    std::vector<supera::LArSpacePoint_t>* _ptr_spacepoint_v;
    std::vector<supera::LArOpFlash_t>* _ptr_opflash_v;
    std::vector<supera::LArCRTHit_t>* _ptr_crthit_v;
    std::string _csv_fname;

    // FIXME(kvtsang) Temporary solution to access associations
    const art::Event* _event;
  };

  //
  // Pointer getter/setters (templated)
  //

  template <>
  const std::vector<supera::LArWire_t>& SuperaBase::LArData<supera::LArWire_t>() const;

  template <>
  const std::vector<supera::LArHit_t>& SuperaBase::LArData<supera::LArHit_t>() const;

  template <>
  const std::vector<supera::LArOpDigit_t>& SuperaBase::LArData<supera::LArOpDigit_t>() const;

  template <>
  const std::vector<supera::LArSimCh_t>& SuperaBase::LArData<supera::LArSimCh_t>() const;

  template <>
  const std::vector<supera::LArMCTruth_t>& SuperaBase::LArData<supera::LArMCTruth_t>() const;

  template <>
  const std::vector<supera::LArMCParticle_t>& SuperaBase::LArData<supera::LArMCParticle_t>() const;

  template <>
  const std::vector<supera::LArMCTrack_t>& SuperaBase::LArData<supera::LArMCTrack_t>() const;

  template <>
  const std::vector<supera::LArMCShower_t>& SuperaBase::LArData<supera::LArMCShower_t>() const;

  template <>
  const std::vector<supera::LArSimEnergyDeposit_t>&
  SuperaBase::LArData<supera::LArSimEnergyDeposit_t>() const;

  template <>
  const std::vector<supera::LArSimEnergyDepositLite_t>&
  SuperaBase::LArData<supera::LArSimEnergyDepositLite_t>() const;

  template <>
  const std::vector<supera::LArSpacePoint_t>& SuperaBase::LArData<supera::LArSpacePoint_t>() const;

  template <>
  const std::vector<supera::LArOpFlash_t>& SuperaBase::LArData<supera::LArOpFlash_t>() const;

  template <>
  const std::vector<supera::LArCRTHit_t>& SuperaBase::LArData<supera::LArCRTHit_t>() const;

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArWire_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArHit_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArOpDigit_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCTruth_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCParticle_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCMiniPart_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCTrack_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArMCShower_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSimCh_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSimEnergyDeposit_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSimEnergyDepositLite_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArSpacePoint_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArOpFlash_t>& data_v);

  template <>
  void SuperaBase::LArData(const std::vector<supera::LArCRTHit_t>& data_v);
  /**
     \class larcv::SuperaBaseFactory
     \brief A concrete factory class for larcv::SuperaBase
  */
  class SuperaBaseProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaBaseProcessFactory() { ProcessFactory::get().add_factory("SuperaBase", this); }
    /// dtor
    ~SuperaBaseProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaBase(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
