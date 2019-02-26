////////////////////////////////////////////////////////////////////////
//
// LArCVMetaMaker.h
//
////////////////////////////////////////////////////////////////////////
#ifndef LARCVMETAMAKER_H
#define LARCVMETAMAKER_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "LArCVMetaData.h"

namespace util{

  class LArCVMetaMaker {

  public:
    LArCVMetaMaker(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    ~LArCVMetaMaker(){};

  public:

    /// Re-configure the service module
    void reconfigure(fhicl::ParameterSet const& pset);

    /// Function to be executed @ run boundary
    void preBeginRun(art::Run const& run);
    void postBeginRun(art::Run const& run);
    /// Function to be executed @ event boundary
    //void preProcessEvent(const art::Event& evt);
    void postProcessEvent(const art::Event& evt);
    /// Function to be executed @ file open
    void preOpenFile(const std::string& filename);
    void postOpenFile(const std::string& filename);

    void postBeginJob();
    void postEndJob();

    std::string GetContent(std::string stream_name) const;

    void addJson(const std::string fname, const std::string strm) { _json_v.push_back(fname + ".json"); _stream_v.push_back(strm); }

  protected:
    bool _quiet;
    std::vector<std::string>              _json_v;
    std::vector<std::string>              _stream_v;
    ::larcv::sam::FCLMetaData_t         _fcl_meta;
    ::larcv::sam::UBMetaData_t          _ub_meta;
    ::larcv::sam::FileCatalogMetaData_t _fcat_meta;
    ::larcv::sam::SAMBuiltInMetaData_t  _sam_meta;

  }; // class LArCVMetaMaker

} //namespace utils

DECLARE_ART_SERVICE(util::LArCVMetaMaker, LEGACY)

#endif 
