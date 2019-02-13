
#include "ubevt/Utilities/FileCatalogMetadataMicroBooNE.h" // RanItay change
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Services/Registry/GlobalSignal.h" // RanItay add
#include "art/Framework/Services/Registry/detail/SignalResponseType.h" // RanItay add
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"



#include "LArCVMetaMaker.h"
#include <TTimeStamp.h>
#include <sstream>
#include <fstream>

//-----------------------------------------------------------------------------------------
util::LArCVMetaMaker::LArCVMetaMaker(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
  : _quiet(true)
//----------------------------------------------------------------------------------------
{
  _json_v.clear();
  reconfigure(pset);

  if(!_quiet) {
    reg.sPostOpenFile.watch     (this, &LArCVMetaMaker::postOpenFile     );
    
    reg.sPreBeginRun.watch      (this, &LArCVMetaMaker::preBeginRun      );
    reg.sPostBeginRun.watch     (this, &LArCVMetaMaker::postBeginRun     );
    
    //reg.sPreProcessEvent.watch  (this, &LArCVMetaMaker::preProcessEvent  );
    //reg.sPostProcessEvent.watch (this, &LArCVMetaMaker::postProcessEvent );
    
    reg.sPostBeginJob.watch     (this, &LArCVMetaMaker::postBeginJob );
    reg.sPostEndJob.watch       (this, &LArCVMetaMaker::postEndJob   );
  }

  _fcl_meta.data_tier = "larcv";
  _fcat_meta.file_format = "larcv";

}

//------------------------------------------------------------------
void util::LArCVMetaMaker::reconfigure(fhicl::ParameterSet const& pset)
//------------------------------------------------------------------
{
  _quiet = !(pset.get<bool>("Enable",false));
}

//------------------------------------------------------------------
void util::LArCVMetaMaker::postBeginJob()
//------------------------------------------------------------------
{
  _sam_meta.start_time = TTimeStamp().AsString("s"); // UTC time stamp in SQL format
  art::ServiceHandle<art::FileCatalogMetadata> larmeta_handle;

  art::FileCatalogMetadata::collection_type larmeta;
  larmeta_handle->getMetadata(larmeta);

  for(auto const& key_value : larmeta) {
    if      (key_value.first == "file_type"         ) _fcat_meta.file_type           = key_value.second;
    else if (key_value.first == "group"             ) _fcat_meta.group               = key_value.second;
    else if (key_value.first == "applicationFamily" ) _fcat_meta.application_family  = key_value.second;
    else if (key_value.first == "process_name"      ) _fcat_meta.application_name    = key_value.second;
    else if (key_value.first == "applicationVersion") _fcat_meta.application_version = key_value.second;
  }
  
}

//------------------------------------------------------------------
void util::LArCVMetaMaker::postEndJob()
//------------------------------------------------------------------
{
  _sam_meta.end_time = TTimeStamp().AsString("s"); // UTC time stamp in SQL format

  art::ServiceHandle<util::FileCatalogMetadataMicroBooNE> ubmeta_handle;

  _fcl_meta.name     = ubmeta_handle->FCLName();
  _fcl_meta.version  = ubmeta_handle->FCLVersion();
  
  _ub_meta.project_name    = ubmeta_handle->ProjectName();
  _ub_meta.project_stage   = ubmeta_handle->ProjectStage();
  _ub_meta.project_version = ubmeta_handle->ProjectVersion();


  //
  // Create json
  //
  for(size_t json_index=0; json_index < _json_v.size(); ++json_index) {

    auto const& json_name = _json_v[json_index];
    auto const& strm_name = _stream_v[json_index];

    std::string content = this->GetContent(strm_name);
    
    std::ofstream fout(json_name);
    fout << content.c_str() << std::endl;
    fout.close();
  }

}

//------------------------------------------------------------
//void util::LArCVMetaMaker::preProcessEvent(const art::Event& evt)
//------------------------------------------------------------
/*{
  std::cout << "This is before event process" << std::endl;
  }*/

//------------------------------------------------------------
void util::LArCVMetaMaker::postProcessEvent(const art::Event& evt)
//------------------------------------------------------------
{
  size_t run    = evt.run();
  size_t subrun = evt.subRun();
  size_t event  = evt.event();

  auto iter = _sam_meta.runs_m.find(run);
  if(iter == _sam_meta.runs_m.end()) {
    // New run found. Insert RunMetaData_t
    iter = (_sam_meta.runs_m.emplace(run,::larcv::sam::RunMetaData_t())).first;
    (*iter).second.run_type = "\" \"";

    art::ServiceHandle<art::FileCatalogMetadata> larmeta_handle;
    art::FileCatalogMetadata::collection_type larmeta;
    larmeta_handle->getMetadata(larmeta);
    for(auto const& key_value : larmeta) {
      if(key_value.first != "run_type") continue;
      (*iter).second.run_type = key_value.second;
      break;
    }
  }
  
  auto& run_meta = (*iter).second;
  run_meta.subruns.insert(subrun);

  if(_sam_meta.event_count == 0) _sam_meta.first_event = event;
  _sam_meta.last_event = event;
  ++_sam_meta.event_count;
}

//------------------------------------------------------
void util::LArCVMetaMaker::preBeginRun(art::Run const& run)
//------------------------------------------------------
{}

//------------------------------------------------------
void util::LArCVMetaMaker::postBeginRun(art::Run const& run)
//------------------------------------------------------
{}

//---------------------------------------------------------------
void util::LArCVMetaMaker::postOpenFile(const std::string& filename)
//---------------------------------------------------------------
{
  _sam_meta.parents.insert(filename);
}

//----------------------------------------------------------------------
std::string util::LArCVMetaMaker::GetContent(std::string stream_name) const
//----------------------------------------------------------------------
{
  std::stringstream msg;
  msg << "{\n"

      << "  \"application\": {\n"
      << "    \"family\"  : " << _fcat_meta.application_family  << ",\n"
      << "    \"name\"    : " << _fcat_meta.application_name    << ",\n"
      << "    \"version\" : " << _fcat_meta.application_version << "\n"
      << "  },\n"
      << "  \"file_format\" : \"" << _fcat_meta.file_format << "\",\n"
      << "  \"file_type\"   : "   << _fcat_meta.file_type   << ",\n"
      << "  \"group\"       : "   << _fcat_meta.group       << ",\n"

      << "  \"fcl.name\"    : \"" << _fcl_meta.name        << "\",\n"
      << "  \"fcl.version\" : \"" << _fcl_meta.version     << "\",\n"
      << "  \"data_tier\"   : \"" << _fcl_meta.data_tier   << "\",\n"
      << "  \"data_stream\" : \"" << stream_name           << "\",\n"

      << "  \"ub_project.name\"    : \"" << _ub_meta.project_name    << "\",\n"
      << "  \"ub_project.stage\"   : \"" << _ub_meta.project_stage   << "\",\n"
      << "  \"ub_project.version\" : \"" << _ub_meta.project_version << "\",\n"

      << "  \"start_time\"  : \"" << _sam_meta.start_time  << "\",\n"
      << "  \"end_time\"    : \"" << _sam_meta.end_time    << "\",\n"
      << "  \"first_event\" : "   << _sam_meta.first_event << ",\n"
      << "  \"last_event\"  : "   << _sam_meta.last_event  << ",\n"
      << "  \"event_count\" : "   << _sam_meta.event_count << ",\n"
      << "  \"parents\"     : [\n";

  size_t ctr=0;
  for(auto const& parent : _sam_meta.parents) {
    size_t name_start = parent.rfind("/");
    if(name_start > parent.length()) name_start = 0;
    else ++name_start;
    msg << "    {  \"file_name\" : \"" << parent.substr(name_start) << "\"  }";
    ++ctr;
    if(ctr==_sam_meta.parents.size()) msg << "\n";
    else msg << ",\n";
  }
  
  msg << "  ],\n"
      << "  \"runs\" : [\n";

  ctr=0;
  for(auto const& run_meta : _sam_meta.runs_m) {
    ++ctr;
    size_t subrun_ctr=0;
    auto const& run = run_meta.first;
    auto const& run_type = run_meta.second.run_type;
    for(auto const& subrun : run_meta.second.subruns) {
      
      msg << "    [  " << run << ",  " << subrun << ",  " << run_type << "]";
      ++subrun_ctr;
      if(ctr == _sam_meta.runs_m.size() && subrun_ctr == run_meta.second.subruns.size())
	msg << "\n";
      else
	msg << ",\n";
    }
  }
  msg << "  ]\n"
      << "}\n";
  return msg.str();
}

namespace util{

  DEFINE_ART_SERVICE(LArCVMetaMaker)

} // namespace util  

