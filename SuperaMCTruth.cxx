#ifndef __SUPERAMCTRUTH_CXX__
#define __SUPERAMCTRUTH_CXX__

#include "SuperaMCTruth.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "Voxel3DSlicer.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "LAr2Image.h"

namespace larcv {

  static SuperaMCTruthProcessFactory __global_SuperaMCTruthProcessFactory__;

  SuperaMCTruth::SuperaMCTruth(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaMCTruth::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsParticle::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    /*
    _store_part = cfg.get<bool>("StoreParticle", true);
    _store_g4_secondary_part = cfg.get<bool>("StoreG4SecondaryParticle", true);
    _store_g4_primary_part = cfg.get<bool>("StoreG4PrimaryParticle", true);
    */
    _pass_origin = cfg.get<unsigned short>("Origin");
    /*
    _filter_pdg  = cfg.get<std::vector<int> >("FilterTargetPDG");
    _filter_min_einit = cfg.get<std::vector<double> >("FilterTargetInitEMin");

    if (_filter_pdg.size() != _filter_min_einit.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetInitEMin not the same length!" << std::endl;
      throw larbys();
    }

    _shower_min_einit = cfg.get<double>("ShowerInitEMin");
    _track_min_einit = cfg.get<double>("TrackInitEMin");
    */
  }

  void SuperaMCTruth::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaMCTruth::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }else if(supera::Voxel3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::Voxel3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());      
    }

    auto& ev_part = mgr.get_data<larcv::EventParticle>(OutParticleLabel());

    auto const& mct_v = LArData<supera::LArMCTruth_t>();
    for(size_t mct_index=0; mct_index<mct_v.size(); ++mct_index) {

      auto const& mct = mct_v[mct_index];

      if(_pass_origin && mct.Origin() != _pass_origin)
	continue;

      for(int i=0; i<mct.NParticles(); ++i) {
	auto const& mcp = mct.GetParticle(i);
	if(mcp.StatusCode() != 1) continue;

	auto const& pos = mcp.Position(0);
	auto const& mom = mcp.Momentum(0);

	/*
	  if ( mcp.PdgCode() != 12 && mcp.PdgCode() != -12 &&
	  mcp.PdgCode() != 14 && mcp.PdgCode() != -14 &&
	  mcp.PdgCode() != 16 && mcp.PdgCode() != -16 )
	  continue;
	*/

	larcv::Particle p;
	p.mct_index(mct_index);
	p.track_id(mcp.TrackId());
	p.shape((mcp.PdgCode() == 11 || mcp.PdgCode() == 22 ? larcv::kShapeShower : larcv::kShapeTrack) );
	p.pdg_code(mcp.PdgCode());
	p.momentum(mom.X()*1.e3,mom.Y()*1.e3,mom.Z()*1.e3);
	p.position(pos.X(),pos.Y(),pos.Z(),pos.T());
	p.energy_init(mom.T());
	p.creation_process(mcp.Process());
	p.parent_track_id(mcp.TrackId());
	p.parent_pdg_code(mcp.PdgCode());
	p.ancestor_track_id(mcp.TrackId());
	p.ancestor_pdg_code(mcp.PdgCode());

	ev_part.emplace_back(std::move(p));
      }
    }
    return true;
  }
      
  void SuperaMCTruth::finalize()
  {}
}

#endif
