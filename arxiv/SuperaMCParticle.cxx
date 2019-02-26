#ifndef __SUPERAMCPARTICLE_CXX__
#define __SUPERAMCPARTICLE_CXX__

#include "SuperaMCParticle.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "Voxel3DSlicer.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "LAr2Image.h"

namespace larcv {

  static SuperaMCParticleProcessFactory __global_SuperaMCParticleProcessFactory__;

  SuperaMCParticle::SuperaMCParticle(const std::string name)
    : SuperaBase(name)
  {}
  /*
  const std::vector<std::pair<size_t, size_t> >& SuperaMCParticle::Particle2MCNode(int part_index) const
  {
    if (part_index >= (int)(_part2mcnode_vv.size()))
      throw larbys("Invalid Particle index requested");

    return _part2mcnode_vv[part_index];
  }
  */
  void SuperaMCParticle::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsParticle::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    /*
    _store_part = cfg.get<bool>("StoreParticle", true);
    _store_g4_secondary_part = cfg.get<bool>("StoreG4SecondaryParticle", true);
    _store_g4_primary_part = cfg.get<bool>("StoreG4PrimaryParticle", true);
    */
    _mcpt.configure(cfg.get<supera::Config_t>("MCParticleTree"));
    _mcpart_maker.configure(cfg.get<supera::Config_t>("MCParticleMaker"));

    _pass_origin = cfg.get<unsigned short>("Origin");
    _mcpt.FilterOrigin(_pass_origin);

    _filter_pdg  = cfg.get<std::vector<int> >("FilterTargetPDG");
    _filter_min_einit = cfg.get<std::vector<double> >("FilterTargetInitEMin");
    _filter_min_edep  = cfg.get<std::vector<double> >("FilterTargetDepEMin");

    if (_filter_pdg.size() != _filter_min_einit.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetInitEMin not the same length!" << std::endl;
      throw larbys();
    }

    if (_filter_pdg.size() != _filter_min_edep.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetDepEMin not the same length!" << std::endl;
      throw larbys();
    }

    _shower_min_einit = cfg.get<double>("ShowerInitEMin");
    _shower_min_edep  = cfg.get<double>("ShowerDepEMin");

    _track_min_einit = cfg.get<double>("TrackInitEMin");
    _track_min_edep  = cfg.get<double>("TrackDepEMin");
    /*
    _filter_min_cols = cfg.get<size_t>("FilterParticleMinCols");
    _filter_min_rows = cfg.get<size_t>("FilterParticleMinRows");
    */
  }

  void SuperaMCParticle::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaMCParticle::process(IOManager& mgr)
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
    bool use_meta3d=Meta3D().valid();

    LARCV_INFO() << "Using " << (use_meta3d ? "3D" : "2D") << " meta..." << std::endl;

    auto ev_part = (EventParticle*)(mgr.get_data("particle", OutParticleLabel()));
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
    //_part2mcnode_vv.clear();

    for (size_t primary_idx = 0; primary_idx < primary_v.size(); ++primary_idx) {

      auto& primary = primary_v[primary_idx];

      // filter out primary of certain origin, if specified
      if (_pass_origin && primary.origin != _pass_origin) {
        LARCV_INFO() << "Skipping a primary " << primary_idx
                     << " (origin type " << primary.origin
                     << " PDG " << primary.pdg
                     << ") with " << primary.daughter_v.size() << " children"
                     << std::endl;
        continue;
      }

      if (!FilterNode(primary)) {
        LARCV_INFO() << "Skipping a primary (TrackID " << primary.track_id
                     << " Origin " << primary.origin
                     << " PDG " << primary.pdg
                     << ") due to FilterNode()" << std::endl;
        continue;
      }
      // at this point we decided to store primary. create one
      _part_v.resize(_part_v.size() + 1);
      auto& pri_part = _part_v.back();

      if(use_meta3d) {
	if (LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	  pri_part = MakeParticle(primary, empty_sch_v, Meta3D());
	else
	  pri_part = MakeParticle(primary, LArData<supera::LArSimCh_t>(), Meta3D());
      }else{
	if (LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	  pri_part = MakeParticle(primary, empty_sch_v);
	else
	  pri_part = MakeParticle(primary, LArData<supera::LArSimCh_t>());
      }
      LARCV_INFO() << "Analyzing primary " << primary_idx << " PDG " << pri_part.pdg_code()
                   << " Origin " << primary.origin
                   << " PDG << " << primary.pdg
                   << " with " << primary.daughter_v.size() << " children" << std::endl;

      std::vector<larcv::Particle> sec_part_v;
      //std::vector<size_t> sec_idx_v;
      for (size_t daughter_idx = 0; daughter_idx < primary.daughter_v.size(); ++daughter_idx) {

        auto const& daughter = primary.daughter_v[daughter_idx];

        if (_pass_origin && daughter.origin != _pass_origin) {
          LARCV_INFO() << "Skipping a daughter " << daughter_idx
                       << " (origin type " << daughter.origin
                       << " PDG " << daughter.pdg
                       << ")"
                       << std::endl;
          continue;
        }

        if (!FilterNode(daughter)) {
          LARCV_INFO() << "Skipping a daughter " << daughter_idx
                       << " (TrackID " << daughter.track_id
                       << " Origin " << daughter.origin
                       << " PDG " << daughter.pdg
                       << ") due to FilterNode()" << std::endl;
          continue;
        }
        larcv::Particle sec_part;
        try {
	  if(use_meta3d) {
	    if (LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	      sec_part = MakeParticle(daughter, empty_sch_v, Meta3D());
	    else
	      sec_part = MakeParticle(daughter, LArData<supera::LArSimCh_t>(), Meta3D());
	  }else{
	    if (LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	      sec_part = MakeParticle(daughter, empty_sch_v);
	    else
	      sec_part = MakeParticle(daughter, LArData<supera::LArSimCh_t>());
	  }
	} catch (const larcv::larbys& err) {
          LARCV_NORMAL() << "Skipping a secondary (PDG,TrackID) = ("
                         << daughter.pdg << "," << daughter.track_id << ") as it could not be turned into Particle" << std::endl;
          continue;
        }
        LARCV_INFO() << "Registering a daughter " << daughter_idx
                     << " (TrackID " << daughter.track_id
                     << " Origin " << daughter.origin
                     << " PDG " << sec_part.pdg_code() << ")" << std::endl;
        sec_part.mct_index(pri_part.mct_index());
        pri_part.energy_deposit(pri_part.energy_deposit() + sec_part.energy_deposit());
        //sec_part_v.emplace_back(std::move(sec_part));
        //sec_idx_v.push_back(daughter_idx);
        _part_v.emplace_back(std::move(sec_part));
      }
      /*
      std::vector<std::pair<size_t, size_t> > part2mcnode_v;
      part2mcnode_v.clear();

      //
      // Register primary Particle
      //
      LARCV_INFO() << "Storing primary Particle (PDG " << primary.part.PdgCode()
                   << " Shape " << primary.part.Shape()
                   << " MCTIndex " << primary.part.mct_index()
                   << " MCSTIndex) " << primary.part.mcst_index()
                   << std::endl;

      _part_v.push_back(primary.part);

      // Record incorporated mcnode
      part2mcnode_v.emplace_back(primary_idx, larcv::kINVALID_SIZE);
      for (size_t daughter_idx = 0; daughter_idx < sec_part_v.size(); ++daughter_idx) {
        auto const& secondary_idx = sec_idx_v[daughter_idx];
        auto const& secondary_part = sec_part_v[daughter_idx];
        if ( (_store_g4_primary_part && secondary_part.TrackID() == secondary_part.ParentTrackID())
             ||
             (_store_g4_secondary_part && secondary_part.TrackID() != secondary_part.ParentTrackID())
           )
          part2mcnode_v.emplace_back(primary_idx, secondary_idx);
      }
      for (size_t daughter_idx = 0; daughter_idx < sec_part_v.size(); ++daughter_idx) {
        auto const& daughter_part = sec_part_v[daughter_idx];
        LARCV_INFO() << "    Associated secondary (PDG " << daughter_part.PdgCode()
                     << " Shape " << daughter_part.Shape()
                     << " MCTIndex " << daughter_part.mct_index()
                     << " MCSTIndex) " << daughter_part.mcst_index()
                     << std::endl;
      }
      _part2mcnode_vv.emplace_back(std::move(part2mcnode_v));

      //
      // Register secondary Particle
      //
      part2mcnode_v.clear();
      for (size_t daughter_idx = 0; daughter_idx < sec_part_v.size(); ++daughter_idx) {
        auto const& secondary_idx = sec_idx_v[daughter_idx];
        auto& secondary_part = sec_part_v[daughter_idx];
        if ( (_store_g4_primary_part && secondary_part.TrackID() == secondary_part.ParentTrackID())
             ||
             (_store_g4_secondary_part && secondary_part.TrackID() != secondary_part.ParentTrackID())
           )
        {
          LARCV_INFO() << "Storing secondary Particle (PDG " << secondary_part.PdgCode()
                       << " Shape " << secondary_part.Shape()
                       << " MCTIndex " << secondary_part.mct_index()
                       << " MCSTIndex) " << secondary_part.mcst_index()
                       << std::endl;
          part2mcnode_v.clear();
          part2mcnode_v.emplace_back(primary_idx, secondary_idx);
          _part_v.emplace_back(std::move(secondary_part));
          _part2mcnode_vv.emplace_back(std::move(part2mcnode_v));
        }
      }
      */
    }
    /*
    if (_store_part)
      ev_part->Emplace(std::move(_part_v));
    */
    ev_part->emplace(std::move(_part_v));
    return true;
  }

  bool SuperaMCParticle::FilterNode(const supera::MCNode & node) const
  {
    if (node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      auto const& mctrack = LArData<supera::LArMCTrack_t>()[node.source_index];
      LARCV_DEBUG() << "MCTrack InitE " << mctrack.Start().E() << " ... DepE "
                    << (mctrack.size() > 1 ? mctrack.front().E() - mctrack.back().E() : 0) << std::endl;
      if (mctrack.Start().E() < _track_min_einit) return false;
      if (_track_min_edep > 0) {
        if (mctrack.size() < 2) return false;
        if ( (mctrack.front().E() - mctrack.back().E()) < _track_min_edep ) return false;
      }
      for (size_t filter_idx = 0; filter_idx < _filter_pdg.size(); ++filter_idx) {
        auto const& pdg = _filter_pdg[filter_idx];
        if (pdg != mctrack.PdgCode()) continue;
        auto const& filter_min_einit = _filter_min_einit[filter_idx];
        if (mctrack.Start().E() < filter_min_einit) return false;
        auto const& filter_min_edep = _filter_min_edep[filter_idx];
        if (filter_min_edep > 0) {
          if (mctrack.size() < 2) return false;
          if ( (mctrack.front().E() - mctrack.back().E()) < filter_min_edep ) return false;
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
        if ( mcshower.DetProfile().E() < filter_min_edep) return false;
      }
    }
    return true;
  }

  larcv::Particle SuperaMCParticle::MakeParticle(const supera::MCNode & node,
						 const std::vector<supera::LArSimCh_t>& sch_v) const
  {
    larcv::Particle res;
    if (node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      auto const& mctrack = LArData<supera::LArMCTrack_t>().at(node.source_index);
      res = _mcpart_maker.MakeParticle(mctrack, sch_v, TimeOffset());
      res.mcst_index(node.source_index);
    }
    else if (node.source_type == supera::MCNode::SourceType_t::kMCShower) {
      auto const& mcshower = LArData<supera::LArMCShower_t>().at(node.source_index);
      res = _mcpart_maker.MakeParticle(mcshower, sch_v, TimeOffset());
      res.mcst_index(node.source_index);
    } else
      throw larbys("Unexpected SourceType_t!");

    // format Particle
    /*
    std::vector<larcv::BBox2D> bb_v;
    for (size_t plane = 0; plane < res.boundingbox_2d().size(); ++plane) {
      auto const& part_meta  = res.boundingbox_2d().at(plane);
      auto const& event_meta = Meta().at(plane);
      bb_v.push_back(FormatMeta(part_meta, event_meta));
    }
    res.boundingbox_2d(bb_v);
    */
    return res;
  }

  larcv::Particle SuperaMCParticle::MakeParticle(const supera::MCNode & node,
						 const std::vector<supera::LArSimCh_t>& sch_v,
						 const Voxel3DMeta& meta3d) const
  {
    larcv::Particle res;
    if (node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      auto const& mctrack = LArData<supera::LArMCTrack_t>().at(node.source_index);
      res = _mcpart_maker.MakeParticle(mctrack, sch_v, TimeOffset(), meta3d);
      res.mcst_index(node.source_index);
    }
    else if (node.source_type == supera::MCNode::SourceType_t::kMCShower) {
      auto const& mcshower = LArData<supera::LArMCShower_t>().at(node.source_index);
      res = _mcpart_maker.MakeParticle(mcshower, sch_v, TimeOffset(), meta3d);
      res.mcst_index(node.source_index);
    } else
      throw larbys("Unexpected SourceType_t!");

    return res;
  }
  /*
  larcv::ImageMeta SuperaMCParticle::FormatMeta(const larcv::ImageMeta & part_image,
      const larcv::ImageMeta & event_image) const
  {

    LARCV_DEBUG() << "Before format  " << part_image.dump();

    const size_t modular_row = event_image.height() / event_image.rows();
    const size_t modular_col = event_image.width()  / event_image.cols();

    if (event_image.rows() < modular_row || event_image.cols() < modular_col) {
      LARCV_ERROR() << "Event image too small to format Particle!" << std::endl;
      throw larbys();
    }
    double min_x  = (part_image.min_x() < event_image.min_x() ? event_image.min_x() : part_image.min_x() );
    double max_y  = (part_image.max_y() > event_image.max_y() ? event_image.max_y() : part_image.max_y() );
    double width  = (part_image.width() + min_x > event_image.max_x() ? event_image.max_x() - min_x : part_image.width());
    double height = (max_y - part_image.height() < event_image.min_y() ? max_y - event_image.min_y() : part_image.height());
    size_t rows   = height / part_image.pixel_height();
    size_t cols   = width  / part_image.pixel_width();

    if (modular_col > 1 && cols % modular_col) {
      int npixels = (modular_col - (cols % modular_col));
      if (event_image.width() < (width + npixels * part_image.pixel_width()))  npixels -= modular_col;
      cols += npixels;
      width += part_image.pixel_width() * npixels;
      if (npixels > 0) {
        // If expanded, make sure it won't go across event_image boundary
        if ( (min_x + width) > event_image.max_x() ) {
          LARCV_DEBUG() << "X: " << min_x << " => " << min_x + width
                        << " exceeds event boundary " << event_image.max_x() << std::endl;
          min_x = event_image.max_x() - width;
        } else if (min_x < event_image.min_x()) {
          LARCV_DEBUG() << "X: " << min_x << " => " << min_x + width
                        << " exceeds event boundary " << event_image.min_x() << std::endl;
          min_x = event_image.min_x();
        }
      }
    }

    if (modular_row > 1 && rows % modular_row) {
      int npixels = (modular_row - (rows % modular_row));
      if (event_image.height() < (height + npixels * part_image.pixel_height()))  npixels -= modular_row;
      rows += npixels;
      height += part_image.pixel_height() * npixels;
      if (npixels > 0) {
        // If expanded, make sure it won't go across event_image boundary
        if ( (max_y - height) < event_image.min_y() ) {
          LARCV_DEBUG() << "Y: " << max_y - height << " => " << max_y
                        << " exceeds event boundary " << event_image.min_y() << std::endl;
          max_y = event_image.min_y() + height;
        } else if (max_y > event_image.max_y()) {
          LARCV_DEBUG() << "Y: " << max_y - height << " => " << max_y
                        << " exceeds event boundary " << event_image.max_y() << std::endl;
          max_y = event_image.max_y();
        }
      }
    }
    LARCV_INFO() << "Creating ImageMeta Width=" << width
                 << " Height=" << height
                 << " NRows=" << (modular_row ? rows / modular_row : rows)
                 << " NCols=" << (modular_col ? cols / modular_col : cols)
                 << " Origin @ (" << min_x << "," << max_y - height << ")" << std::endl;
    larcv::ImageMeta res(min_x, max_y - height, min_x + width, max_y,
                         cols, rows,
                         part_image.id(), larcv::kUnitWireTime);

    LARCV_DEBUG() << "Event image   " << event_image.dump();

    LARCV_DEBUG() << "After format  " << res.dump();

    //res = event_image.overlap(res);
    //LARCV_DEBUG() << "After overlap " << res.dump();
    return res;
  }
  */
  void SuperaMCParticle::finalize()
  {}
}

#endif
