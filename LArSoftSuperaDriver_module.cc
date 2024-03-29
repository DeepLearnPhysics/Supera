////////////////////////////////////////////////////////////////////////
// Class:       LArSoftSuperaDriver
// Module Type: analyzer
// File:        LArSoftSuperaDriver_module.cc
//
// Generated at Tue May  2 17:36:15 2017 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "FMWKInterface.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include <TString.h>
#include <TTimeStamp.h>

//#include "LArCVMetaMaker.h"
#include "LArCVSuperaDriver.h"
#include "GenRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "TRandom.h"
#include "nurandom/RandomUtils/NuRandomService.h"

class LArSoftSuperaDriver;

class LArSoftSuperaDriver : public art::EDAnalyzer {
public:
  explicit LArSoftSuperaDriver(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LArSoftSuperaDriver(LArSoftSuperaDriver const &) = delete;
  LArSoftSuperaDriver(LArSoftSuperaDriver &&) = delete;
  LArSoftSuperaDriver & operator = (LArSoftSuperaDriver const &) = delete;
  LArSoftSuperaDriver & operator = (LArSoftSuperaDriver &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob() override;

  template<class LArSoftDataType> void get_label(const art::Event& e, ::supera::LArDataType_t SuperaDataType, bool checkLength = false);

private:

  // Declare member data here.
  larcv::LArCVSuperaDriver _supera;
  unsigned int _verbosity;
  CLHEP::HepRandomEngine& fFlatEngine;
  bool _strictDataLoading;
};


LArSoftSuperaDriver::LArSoftSuperaDriver(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
  , fFlatEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
                  createEngine(0, "HepJamesRandom", "Supera"), "HepJamesRandom", "Supera"))
 // More initializers here.
{
  //fEfficiencyEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "Efficiencies"))
  _verbosity = p.get<unsigned int>("Verbosity",3);
  _strictDataLoading = p.get<bool>("StrictDataLoading", true);
  // Setup random number generator
  supera::GenRandom::get().SetFlatGen(new CLHEP::RandFlat(fFlatEngine,0,1));

  std::string supera_cfg;
  cet::search_path finder("FHICL_FILE_PATH");

  if( !finder.find_file(p.get<std::string>("supera_params"), supera_cfg) )
    throw cet::exception("LArSoftSuperaDriver") << "Unable to find supera cfg in "  << finder.to_string() << "\n";

  _supera.configure(supera_cfg);

  _supera.set_verbosity(::larcv::msg::Level_t(_verbosity));

  auto process_names = _supera.ProcessNames();
  for(auto const& proc_name : process_names) {
    std::string param_name = "CSV" + proc_name;
    auto constraint_file = p.get<std::string>(param_name.c_str(),"");
    if(constraint_file.empty()) continue;

    std::string fullpath;
    if( !finder.find_file(constraint_file,fullpath) )
      throw cet::exception("LArSoftSuperaDriver") << "Unable to find CSV file "  << constraint_file << "\n";

    _supera.SetCSV(proc_name,fullpath);
  }

  // Decide on output filename
  auto out_fname = p.get<std::string>("out_filename","");
  if(p.get<bool>("unique_filename") || out_fname.empty()) {
    TString tmp_fname;
    if(out_fname.empty())
      tmp_fname = "emptyname_" + p.get<std::string>("stream");
    else {
      tmp_fname = out_fname;
      tmp_fname.ReplaceAll(".root","");
      tmp_fname += "_" + p.get<std::string>("stream");
    }
    TTimeStamp ts;
    out_fname = Form("%s_%08d_%06d_%06d.root",tmp_fname.Data(),ts.GetDate(),ts.GetTime(), (int)(ts.GetNanoSec()/1.e3));
  }

  _supera.override_output_file(out_fname);

  //art::ServiceHandle<util::LArCVMetaMaker> metamaker;
  //metamaker->addJson(out_fname,p.get<std::string>("stream"));
}

void LArSoftSuperaDriver::beginJob()
{
  _supera.initialize();
}

// Define boilerplate function to be used for
// various data types. Use templates for now.
template <class LArSoftDataType> void LArSoftSuperaDriver::get_label(const art::Event& e, ::supera::LArDataType_t SuperaDataType, bool checkLength) {
  for(auto const& label : _supera.DataLabels(SuperaDataType)) {
    if(label.empty()) continue;
    art::Handle<std::vector<LArSoftDataType> > data_h;
    if(label.find(" ")<label.size()) {
      e.getByLabel(label.substr(0,label.find(" ")),
       label.substr(label.find(" ")+1,label.size()-label.find(" ")-1),
       data_h);
    }else{ e.getByLabel(label, data_h); }
    if(!data_h.isValid() || (checkLength ? data_h->empty() : false)) {
      std::cerr<< "Attempted to load data: " << label << std::endl;
      if(_strictDataLoading)
        throw ::larcv::larbys("Could not locate data!");
      else
        return;
    }
    _supera.SetDataPointer(*data_h,label);
  }
}

void LArSoftSuperaDriver::analyze(art::Event const & e)
{
  // FIXME(kvtsang) Temporary solution to access associations
  _supera.SetEvent(&e);

  //
  // set data pointers
  //



  // hit
  //get_label<::supera::LArDataType_t::kLArHit_t, recob::Hit>();

  // wire
  if(_verbosity==0) std::cout << "Checking Wire data request" << std::endl;
  get_label<recob::Wire>(e, ::supera::LArDataType_t::kLArWire_t);

  // opdigit
  //get_label<::supera::LArDataType_t::kLArOpDigit_t, raw::OpDetWaveform>();

  // mctruth
  if(_verbosity==0) std::cout << "Checking MCTruth data request" << std::endl;
  get_label<simb::MCTruth>(e, ::supera::LArDataType_t::kLArMCTruth_t);

  // mcparticle
  if(_verbosity==0) std::cout << "Checking MCParticle data request" << std::endl;
  get_label<simb::MCParticle>(e, ::supera::LArDataType_t::kLArMCParticle_t, true);

  // mcminipart
  if(_verbosity==0) std::cout << "Checking MCMiniPart data request" << std::endl;
  get_label<sim::MCParticleLite>(e, ::supera::LArDataType_t::kLArMCMiniPart_t, true);
  //
  // SimEnergyDeposit
  if(_verbosity==0) std::cout << "Checking SimEnergyDeposit data request" << std::endl;
  get_label<sim::SimEnergyDeposit>(e, ::supera::LArDataType_t::kLArSimEnergyDeposit_t);

  // SimEnergyDepositLite
  if(_verbosity==0) std::cout << "Checking SimEnergyDepositLite data request" << std::endl;
  get_label<sim::SimEnergyDepositLite>(e, ::supera::LArDataType_t::kLArSimEnergyDepositLite_t);

  // mctrack
  get_label<sim::MCTrack>(e, ::supera::LArDataType_t::kLArMCTrack_t);

  // mcshower
  get_label<sim::MCShower>(e, ::supera::LArDataType_t::kLArMCShower_t);

  // SpacePoint
  get_label<recob::SpacePoint>(e, ::supera::LArDataType_t::kLArSpacePoint_t);

  // simch
  get_label<sim::SimChannel>(e, ::supera::LArDataType_t::kLArSimCh_t);

  // OpFlash
  if(_verbosity==0) std::cout << "Checking OpFlash data request" << std::endl;
  get_label<recob::OpFlash>(e, ::supera::LArDataType_t::kLArOpFlash_t);

  // CRTHit
  if(_verbosity==0) std::cout << "Checking OpFlash data request" << std::endl;
  get_label<sbn::crt::CRTHit>(e, ::supera::LArDataType_t::kLArCRTHit_t);

  /*
  // chstatus
  auto supera_chstatus = _supera.SuperaChStatusPointer();
  if(supera_chstatus) {

    // Set database status
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    for(size_t i=0; i < geom->Nchannels(); ++i) {
      auto const wid = geom->ChannelToWire(i).front();
      if (!chanFilt.IsPresent(i)) supera_chstatus->set_chstatus(wid.Plane, wid.Wire, ::larcv::chstatus::kNOTPRESENT);
      else supera_chstatus->set_chstatus(wid.Plane, wid.Wire, (short)(chanFilt.Status(i)));
    }
  }
  */

  //
  // execute supera
  //
  _supera.process(e.id().run(),e.id().subRun(),e.id().event());
}

void LArSoftSuperaDriver::endJob()
{
  _supera.finalize();
}

DEFINE_ART_MODULE(LArSoftSuperaDriver)
