#ifndef __SUPERAMCTRUTH_CXX__
#define __SUPERAMCTRUTH_CXX__

#include "SuperaMCTruth.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventNeutrino.h"
#include "larcv/core/DataFormat/Neutrino.h"

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
    _producer_labels = cfg.get<std::vector<std::string>>("MCTruthProducers", {});

    // Backward compatibility
    if (_producer_labels.empty()) {
      auto prod = cfg.get<std::string>("LArMCTruthProducer");
      _producer_labels.push_back(prod);
    }

    for (auto const& label : _producer_labels)
        Request(supera::LArDataType_t::kLArMCTruth_t, label);
  }

  void SuperaMCTruth::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaMCTruth::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto& ev_part = mgr.get_data<larcv::EventParticle>(_output_label);
		auto& ev_nu   = mgr.get_data<larcv::EventNeutrino>(_output_label);

    auto const *ev = GetEvent();
    size_t label_idx = 0;
    for (auto const& label : _producer_labels) {
      //std::cout << "looping for " << label << std::endl;
      auto handle = ev->getHandle<std::vector<simb::MCTruth>>(label);

      if (! handle.isValid()) {
        std::cerr << "Failed to get MCTruth from " << label << std::endl;
          //return false;
          continue;
      }
      //auto const& mct_v = LArData<supera::LArMCTruth_t>();
      auto mct_v = *handle;
      for(size_t mct_index=0; mct_index<mct_v.size(); ++mct_index) {
      //std::cout << "mct_index is " << label << " " << label_idx << std::endl;

      auto const& mct = mct_v[mct_index];

      if(_pass_origin && mct.Origin() != _pass_origin)
	continue;

      if(mct.NeutrinoSet()) {
	auto const& mcnu = mct.GetNeutrino().Nu();
	auto const& mcnuint = mct.GetNeutrino();

	larcv::Particle nu;
	//nu.mct_index(mct_index); // this is always 0 ? when would there be several MCTruth for same generator?
	nu.mct_index(label_idx); // for now store which generator this came from
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
	nu.nu_interaction_type(mcnuint.InteractionType());
	nu.nu_current_type(mcnuint.CCNC());
	//nu.nu_interaction_mode(mcnuint.Mode());
	//nu.nu_nucleon(mcnuint.HitNuc());
	//nu.nu_quark(mcnuint.HitQuark());
	//nu.nu_momentum_transfer(mcnuint.QSqr());
	//nu.nu_bjorken_x(mcnuint.X());

	LARCV_INFO() << nu.dump() << std::endl;
	LARCV_INFO() << mcnuint << std::endl;
	ev_part.emplace_back(std::move(nu));

	larcv::Neutrino neutrino;
  //neutrino.mct_index(mct_index);
	neutrino.mct_index(label_idx);
	neutrino.nu_track_id(mcnu.TrackId());
	neutrino.lepton_track_id(mcnuint.Lepton().TrackId());
	neutrino.current_type(mcnuint.CCNC());
	neutrino.interaction_mode(mcnuint.Mode());
	neutrino.interaction_type(mcnuint.InteractionType());
	neutrino.target(mcnuint.Target());
	neutrino.nucleon(mcnuint.HitNuc());
	neutrino.quark(mcnuint.HitQuark());
	neutrino.hadronic_invariant_mass(mcnuint.W());
	neutrino.bjorken_x(mcnuint.X());
	neutrino.inelasticity(mcnuint.Y());
	neutrino.momentum_transfer(mcnuint.QSqr());
	neutrino.theta(mcnuint.Theta());
	neutrino.pdg_code(mcnu.PdgCode());
	neutrino.momentum(mcnu.Px(), mcnu.Py(), mcnu.Pz());
	neutrino.position(mcnu.Vx(), mcnu.Vy(), mcnu.Vz(), mcnu.T());
	neutrino.energy_init(mcnu.Momentum(0).T());
	neutrino.creation_process(mcnu.Process());
	auto& traj = mcnu.Trajectory();
	//auto& processes = traj.TrajectoryProcesses();
	for(size_t j = 0; j < traj.size(); ++j) {
		neutrino.add_trajectory_point(
			traj.X(j), traj.Y(j), traj.Z(j), traj.T(j),
			traj.Px(j), traj.Py(j), traj.Pz(j), traj.E(j)
		);
	}
	ev_nu.emplace_back(std::move(neutrino));
      }

      for(int i=0; i<mct.NParticles(); ++i) {
	auto const& mcp = mct.GetParticle(i);
	if(mcp.StatusCode() != 1) std::cout << "bad particle: track id = " << mcp.TrackId() << " PDG = " << mcp.PdgCode() << " StatusCode = " << mcp.StatusCode() << std::endl;
	if(mcp.StatusCode() != 1) continue;

	auto const& pos = mcp.Position(0);
	auto const& mom = mcp.Momentum(0);

	larcv::Particle p;
  //p.mct_index(mct_index);
	p.mct_index(label_idx);
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
	LARCV_DEBUG() << "*******------ track id " << mcp.TrackId() << std::endl;
	LARCV_DEBUG() << "*******------ ancestor pdg " << mcp.PdgCode() << std::endl;
	LARCV_DEBUG() << "*******------ ancestor id " << mcp.TrackId() << std::endl;

	LARCV_INFO() << p.dump() << std::endl;

	ev_part.emplace_back(std::move(p));
      } // end of particle loop
    } // end of mct loop
    label_idx++;
    } // end of producer loop
    return true;
  }

  void SuperaMCTruth::finalize()
  {}
}

#endif
