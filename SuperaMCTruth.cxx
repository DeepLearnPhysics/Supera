#ifndef __SUPERAMCTRUTH_CXX__
#define __SUPERAMCTRUTH_CXX__

#include "SuperaMCTruth.h"
#include "larcv/core/DataFormat/EventParticle.h"

namespace larcv {

  static SuperaMCTruthProcessFactory __global_SuperaMCTruthProcessFactory__;

  SuperaMCTruth::SuperaMCTruth(const std::string name)
    : SuperaBase(name)
  {}
  void SuperaMCTruth::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _output_label = cfg.get<std::string>("OutParticleLabel");
    _pass_origin = cfg.get<unsigned short>("Origin");
  }

  void SuperaMCTruth::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaMCTruth::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto& ev_part = mgr.get_data<larcv::EventParticle>(_output_label);

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
