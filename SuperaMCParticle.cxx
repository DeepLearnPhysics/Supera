#ifndef __SUPERAMCPARTICLE_CXX__
#define __SUPERAMCPARTICLE_CXX__

#include "SuperaMCParticle.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaMCParticleProcessFactory __global_SuperaMCParticleProcessFactory__;

  SuperaMCParticle::SuperaMCParticle(const std::string name) : SuperaBase(name) {}

  void SuperaMCParticle::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);

    _output_label = cfg.get<std::string>("OutParticleLabel");
    _ref_meta3d_cluster3d = cfg.get<std::string>("Meta3DFromCluster3D", "");
    _ref_meta3d_tensor3d = cfg.get<std::string>("Meta3DFromTensor3D", "");

    _mcpt.configure(cfg.get<supera::Config_t>("MCParticleTree"));
    _mcpart_maker.configure(cfg.get<supera::Config_t>("MCParticleMaker"));

    _pass_origin = cfg.get<unsigned short>("Origin");
    _mcpt.FilterOrigin(_pass_origin);

    _filter_pdg = cfg.get<std::vector<int>>("FilterTargetPDG");
    _filter_min_einit = cfg.get<std::vector<double>>("FilterTargetInitEMin");
    _filter_min_edep = cfg.get<std::vector<double>>("FilterTargetDepEMin");

    if (_filter_pdg.size() != _filter_min_einit.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetInitEMin not the same length!"
                       << std::endl;
      throw larbys();
    }

    if (_filter_pdg.size() != _filter_min_edep.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetDepEMin not the same length!"
                       << std::endl;
      throw larbys();
    }

    _shower_min_einit = cfg.get<double>("ShowerInitEMin");
    _shower_min_edep = cfg.get<double>("ShowerDepEMin");

    _track_min_einit = cfg.get<double>("TrackInitEMin");
    _track_min_edep = cfg.get<double>("TrackDepEMin");
  }

  void SuperaMCParticle::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaMCParticle::process(IOManager& mgr)
  {
    LARCV_DEBUG() << "****---- process MCParticle" << std::endl;
    SuperaBase::process(mgr);

    auto const& part_v = LArData<supera::LArMCParticle_t>();
    std::vector<int> track_ids;
    track_ids.reserve(part_v.size());
    for (auto const& p : part_v) {
      if (abs(p.TrackId()) >= track_ids.size()) track_ids.resize(abs(p.TrackId()) + 1);
      track_ids[p.TrackId()] = abs(p.PdgCode());
    }
    for (auto const& p : part_v) {
      if (abs(p.PdgCode()) != 11) continue;
      if (p.Process() != "muIoni") continue;
      if (track_ids[p.Mother()] != 13) continue;
      LARCV_DEBUG() << "LOGME " << p.E() * 1000. << std::endl;
    }

    return true;
    larcv::Voxel3DMeta meta;
    if (!_ref_meta3d_cluster3d.empty()) {
      auto const& ev_cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(_ref_meta3d_cluster3d);
      meta = ev_cluster3d.meta();
    }
    else if (!_ref_meta3d_tensor3d.empty()) {
      auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_ref_meta3d_tensor3d);
      meta = ev_tensor3d.meta();
    }

    auto ev_part = (EventParticle*)(mgr.get_data("particle", _output_label));
    if (!ev_part) {
      LARCV_CRITICAL() << "Output part could not be created!" << std::endl;
      throw larbys();
    }
    if (!(ev_part->as_vector().empty())) {
      LARCV_CRITICAL() << "Output part array not empty!" << std::endl;
      throw larbys();
    }

    LARCV_INFO() << "Running MCParticleTree::Register" << std::endl;
    _mcpt.Register(LArData<supera::LArMCTrack_t>(), LArData<supera::LArMCShower_t>());
    //_mcpt.dump();

    auto primary_v = _mcpt.PrimaryArray();
    LARCV_INFO() << "Found " << primary_v.size() << " primary particles" << std::endl;

    std::vector<supera::LArSimCh_t> empty_sch_v;

    _part_v.clear();

    for (size_t primary_idx = 0; primary_idx < primary_v.size(); ++primary_idx) {

      auto& primary = primary_v[primary_idx];

      // filter out primary of certain origin, if specified
      if (_pass_origin && primary.origin != _pass_origin) {
        LARCV_INFO() << "Skipping a primary " << primary_idx << " (origin type " << primary.origin
                     << " PDG " << primary.pdg << ") with " << primary.daughter_v.size()
                     << " children" << std::endl;
        continue;
      }

      if (!FilterNode(primary)) {
        LARCV_INFO() << "Skipping a primary (TrackID " << primary.track_id << " Origin "
                     << primary.origin << " PDG " << primary.pdg << ") due to FilterNode()"
                     << std::endl;
        continue;
      }
      // at this point we decided to store primary. create one
      //_part_v.resize(_part_v.size() + 1);
      //auto& pri_part = _part_v.back();
      LARCV_DEBUG() << "****---- pri_part" << std::endl;
      auto pri_part = MakeParticle(primary, meta);
      _part_v.push_back(std::move(pri_part));

      LARCV_INFO() << "Analyzing primary " << primary_idx << " PDG " << pri_part.pdg_code()
                   << " Origin " << primary.origin << " with " << primary.daughter_v.size()
                   << " children" << std::endl;

      std::vector<larcv::Particle> sec_part_v;
      for (size_t daughter_idx = 0; daughter_idx < primary.daughter_v.size(); ++daughter_idx) {

        auto const& daughter = primary.daughter_v[daughter_idx];

        if (_pass_origin && daughter.origin != _pass_origin) {
          LARCV_INFO() << "Skipping a daughter " << daughter_idx << " (origin type "
                       << daughter.origin << " PDG " << daughter.pdg << ")" << std::endl;
          continue;
        }

        if (!FilterNode(daughter)) {
          LARCV_INFO() << "Skipping a daughter " << daughter_idx << " (TrackID "
                       << daughter.track_id << " Origin " << daughter.origin << " PDG "
                       << daughter.pdg << ") due to FilterNode()" << std::endl;
          continue;
        }
        larcv::Particle sec_part;
        try {
          LARCV_DEBUG() << "****---- sec_part" << std::endl;
          sec_part = MakeParticle(daughter, meta);
        }
        catch (const larcv::larbys& err) {
          LARCV_NORMAL() << "Skipping a secondary (PDG,TrackID) = (" << daughter.pdg << ","
                         << daughter.track_id << ") as it could not be turned into Particle"
                         << std::endl;
          continue;
        }
        LARCV_INFO() << "Registering a daughter " << daughter_idx << " (TrackID "
                     << daughter.track_id << " Origin " << daughter.origin << " PDG "
                     << sec_part.pdg_code() << ")" << std::endl;
        sec_part.mct_index(pri_part.mct_index());
        pri_part.energy_deposit(pri_part.energy_deposit() + sec_part.energy_deposit());
        _part_v.push_back(std::move(sec_part));
      }
    }

    ev_part->emplace(std::move(_part_v));
    return true;
  }

  bool SuperaMCParticle::FilterNode(const supera::MCNode& node) const
  {
    if (node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      auto const& mctrack = LArData<supera::LArMCTrack_t>()[node.source_index];
      LARCV_DEBUG() << "MCTrack InitE " << mctrack.Start().E() << " ... DepE "
                    << (mctrack.size() > 1 ? mctrack.front().E() - mctrack.back().E() : 0)
                    << std::endl;
      if (mctrack.Start().E() < _track_min_einit) return false;
      if (_track_min_edep > 0) {
        if (mctrack.size() < 2) return false;
        if ((mctrack.front().E() - mctrack.back().E()) < _track_min_edep) return false;
      }
      for (size_t filter_idx = 0; filter_idx < _filter_pdg.size(); ++filter_idx) {
        auto const& pdg = _filter_pdg[filter_idx];
        if (pdg != mctrack.PdgCode()) continue;
        auto const& filter_min_einit = _filter_min_einit[filter_idx];
        if (mctrack.Start().E() < filter_min_einit) return false;
        auto const& filter_min_edep = _filter_min_edep[filter_idx];
        if (filter_min_edep > 0) {
          if (mctrack.size() < 2) return false;
          if ((mctrack.front().E() - mctrack.back().E()) < filter_min_edep) return false;
        }
      }
    }
    else if (node.source_type == supera::MCNode::SourceType_t::kMCShower) {
      auto const& mcshower = LArData<supera::LArMCShower_t>()[node.source_index];
      LARCV_DEBUG() << "MCShower InitE " << mcshower.Start().E() << " ... DepE "
                    << mcshower.DetProfile().E() << std::endl;
      if (mcshower.Start().E() < _shower_min_einit) return false;
      if (mcshower.DetProfile().E() < _shower_min_edep) return false;
      for (size_t filter_idx = 0; filter_idx < _filter_pdg.size(); ++filter_idx) {
        auto const& pdg = _filter_pdg[filter_idx];
        if (pdg != mcshower.PdgCode()) continue;
        auto const& filter_min_einit = _filter_min_einit[filter_idx];
        if (mcshower.Start().E() < filter_min_einit) return false;
        auto const& filter_min_edep = _filter_min_edep[filter_idx];
        if (mcshower.DetProfile().E() < filter_min_edep) return false;
      }
    }
    return true;
  }

  larcv::Particle SuperaMCParticle::MakeParticle(const supera::MCNode& node,
                                                 const larcv::Voxel3DMeta& meta) const
  {
    LARCV_DEBUG() << "***--- SuperaMCParticle::MakeParticle" << std::endl;
    larcv::Particle res;
    if (node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      LARCV_DEBUG() << "**-- supera::MCNode::SourceType_t::kMCTrack" << std::endl;
      auto const& mctrack = LArData<supera::LArMCTrack_t>().at(node.source_index);
      res = _mcpart_maker.MakeParticle(mctrack, meta);
      res.mcst_index(node.source_index);
    }
    else if (node.source_type == supera::MCNode::SourceType_t::kMCShower) {
      LARCV_DEBUG() << "**-- supera::MCNode::SourceType_t::kMCTrack" << std::endl;
      auto const& mcshower = LArData<supera::LArMCShower_t>().at(node.source_index);
      res = _mcpart_maker.MakeParticle(mcshower);
      res.mcst_index(node.source_index);
    }
    else
      throw larbys("Unexpected SourceType_t!");

    return res;
  }

  void SuperaMCParticle::finalize() {}
}

#endif
