#ifndef __SUPERAMETAMAKER_CXX__
#define __SUPERAMETAMAKER_CXX__

#include "SuperaMetaMaker.h"
#include "SuperaCSVReader.h"
#include "LAr2Image.h"
#include "ImageMetaMakerFactory.h"
#include "ImageMetaMaker.h"
#include "PulledPork3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"

namespace larcv {

  static SuperaMetaMakerProcessFactory __global_SuperaMetaMakerProcessFactory__;

  SuperaMetaMaker::SuperaMetaMaker(const std::string name)
    : SuperaBase(name)
    , _meta_maker(nullptr)
  {}
  
  bool SuperaMetaMaker::is(const std::string question) const
  {
    if(question == "Supera") return true;
    if(question == "SuperaMetaMaker") return true;
    return false;
  }

  void SuperaMetaMaker::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    if(_meta_maker) delete _meta_maker;
    _meta_maker = supera::CreateImageMetaMaker(cfg);
    supera::ImageMetaMaker::SetSharedMetaMaker(_meta_maker);
  }

  void SuperaMetaMaker::initialize()
  { SuperaBase::initialize(); }

  bool SuperaMetaMaker::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    if(supera::PulledPork3DSlicer::Is(_meta_maker)) {
      ((supera::PulledPork3DSlicer*)(_meta_maker))->ClearEventData();
      if(!LArDataLabel(supera::LArDataType_t::kLArMCTruth_t).empty())
	((supera::PulledPork3DSlicer*)(_meta_maker))->AddConstraint(LArData<supera::LArMCTruth_t>());
      if(!CSV().empty() && _constraint_m.empty()) {
	_constraint_m.clear();
	supera::csvreader::read_constraint_file(CSV(),_constraint_m);
	LARCV_NORMAL() << "Loaded constraint points for " << _constraint_m.size() << " events..." << std::endl;
      }
      if(_constraint_m.size()) {
	auto const& larcv_event_id = mgr.event_id();
	//supera::RSEID supera_event_id(larcv_event_id.run(),larcv_event_id.subrun(),larcv_event_id.event());
	supera::RSEID supera_event_id(larcv_event_id.run(),larcv_event_id.event());
	auto iter = _constraint_m.find(supera_event_id);
	if(iter!=_constraint_m.end()) 
	  ((supera::PulledPork3DSlicer*)(_meta_maker))->AddConstraint((*iter).second[0],(*iter).second[1],(*iter).second[2]);
      }
      if(!LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	((supera::PulledPork3DSlicer*)(_meta_maker))->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
      else
	((supera::PulledPork3DSlicer*)(_meta_maker))->GenerateMeta(TimeOffset());
    }
    
    return true;
  }

  void SuperaMetaMaker::finalize()
  {}


}
#endif
