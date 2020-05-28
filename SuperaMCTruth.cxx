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

      if(mct.NeutrinoSet()) {
	auto const& mcnu = mct.GetNeutrino().Nu();
	larcv::Particle nu;
	nu.mct_index(mct_index);
	nu.track_id(mcnu.TrackId());
	nu.pdg_code(mcnu.PdgCode());
	nu.momentum(mcnu.Momentum(0).X()*1.e3,
		    mcnu.Momentum(0).Y()*1.e3,
		    mcnu.Momentum(0).Z()*1.e3);
	nu.position(mcnu.Position(0).X(),
		    mcnu.Position(0).Y(),
		    mcnu.Position(0).Z(),
		    mcnu.Position(0).T());
	nu.energy_init(mcnu.Momentum(0).T());
	nu.creation_process(mcnu.Process());
	nu.parent_track_id(mcnu.TrackId());
	nu.parent_pdg_code(mcnu.PdgCode());
	nu.ancestor_track_id(mcnu.TrackId());
	nu.ancestor_pdg_code(mcnu.PdgCode());
	
	LARCV_INFO() << nu.dump() << std::endl;
	ev_part.emplace_back(std::move(nu));
      }

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

	LARCV_INFO() << p.dump() << std::endl;

	ev_part.emplace_back(std::move(p));
      }
    }
    return true;
  }
      
  void SuperaMCTruth::finalize()
  {}
}

#endif
