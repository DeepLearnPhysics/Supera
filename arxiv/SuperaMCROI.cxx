#ifndef __SUPERAMCPARTICLE_CXX__
#define __SUPERAMCPARTICLE_CXX__

#include "SuperaMCParticle.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "LAr2Image.h"

namespace larcv {

  static SuperaMCParticleProcessFactory __global_SuperaMCParticleProcessFactory__;

  SuperaMCParticle::SuperaMCParticle(const std::string name)
    : SuperaBase(name)
  {}

  const std::vector<std::pair<size_t,size_t> >& SuperaMCParticle::Particle2MCNode(int part_index) const
  {
    if(part_index >= (int)(_part2mcnode_vv.size()))
      throw larbys("Invalid Particle index requested");

    return _part2mcnode_vv[part_index];
  }
    
  void SuperaMCParticle::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsParticle::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    
    _store_part = cfg.get<bool>("StoreParticle",true);
    _store_g4_secondary_part = cfg.get<bool>("StoreG4SecondaryParticle",true);
    _store_g4_primary_part = cfg.get<bool>("StoreG4PrimaryParticle",true);
    _mcpt.configure(cfg.get<supera::Config_t>("MCParticleTree"));
    _mcpart_maker.configure(cfg.get<supera::Config_t>("MCParticleMaker"));

    _pass_origin = cfg.get<unsigned short>("Origin");
    _mcpt.FilterOrigin(_pass_origin);
    
    _filter_pdg  = cfg.get<std::vector<int> >("FilterTargetPDG");
    _filter_min_einit = cfg.get<std::vector<double> >("FilterTargetInitEMin");
    _filter_min_edep  = cfg.get<std::vector<double> >("FilterTargetDepEMin");

    if(_filter_pdg.size() != _filter_min_einit.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetInitEMin not the same length!" << std::endl;
      throw larbys();
    }

    if(_filter_pdg.size() != _filter_min_edep.size()) {
      LARCV_CRITICAL() << "FilterTargetPDG and FilterTargetDepEMin not the same length!" << std::endl;
      throw larbys();
    }

    _shower_min_einit = cfg.get<double>("ShowerInitEMin");
    _shower_min_edep  = cfg.get<double>("ShowerDepEMin");

    _track_min_einit = cfg.get<double>("TrackInitEMin");
    _track_min_edep  = cfg.get<double>("TrackDepEMin");

    _filter_min_cols = cfg.get<size_t>("FilterParticleMinCols");
    _filter_min_rows = cfg.get<size_t>("FilterParticleMinRows");
    
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
    }
    
    auto const& meta_v = Meta();

    if(meta_v.empty()) {
      LARCV_CRITICAL() << "Meta not created!" << std::endl;
      throw larbys();
    }
    auto ev_part = (EventParticle*)(mgr.get_data("part",OutParticleLabel()));
    if(!ev_part) {
      LARCV_CRITICAL() << "Output part could not be created!" << std::endl;
      throw larbys();
    }
    if(!(ev_part->ParticleArray().empty())) {
      LARCV_CRITICAL() << "Output part array not empty!" << std::endl;
      throw larbys();
    }

    LARCV_INFO() << "Running MCParticleTree::Register" << std::endl;
    _mcpt.Register(LArData<supera::LArMCTrack_t>(),LArData<supera::LArMCShower_t>());
    _mcpt.dump();
    
    auto primary_v = _mcpt.PrimaryArray();
    LARCV_INFO() << "Found " << primary_v.size() << " primary particles" << std::endl;

    std::vector<supera::LArSimCh_t> empty_sch_v;

    _part_v.clear();
    _part2mcnode_vv.clear();

    for(size_t primary_idx = 0; primary_idx < primary_v.size(); ++primary_idx) {
      
      auto& primary = primary_v[primary_idx];

      // filter out primary of certain origin, if specified
      if(_pass_origin && primary.origin != _pass_origin) {
	LARCV_INFO() << "Skipping a primary " << primary_idx
		     << " (origin type " << primary.origin
		     << " PDG " << primary.pdg
		     << ") with " << primary.daughter_v.size() << " children"
		     << std::endl;
	continue;
      }

      if(!FilterNode(primary)) {
	LARCV_INFO() << "Skipping a primary (TrackID " << primary.track_id
		     << " Origin " << primary.origin
		     << " PDG " << primary.pdg
		     << ") due to FilterNode()" << std::endl;
	continue;
      }
      
      if(LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	primary.part = MakeParticle(primary,empty_sch_v);	  
      else
	primary.part = MakeParticle(primary,LArData<supera::LArSimCh_t>());

      LARCV_INFO() << "Analyzing primary " << primary_idx << " PDG " << primary.part.PdgCode()
		   << " Origin " << primary.origin
		   << " PDG << " << primary.pdg
		   << " with " << primary.daughter_v.size() << " children" << std::endl;

      std::vector<larcv::Particle> sec_part_v;
      std::vector<size_t> sec_idx_v;
      for(size_t daughter_idx=0; daughter_idx<primary.daughter_v.size(); ++daughter_idx) {

	auto const& daughter = primary.daughter_v[daughter_idx];

	if(!FilterNode(daughter)) {
	  LARCV_INFO() << "Skipping a daughter " << daughter_idx
		       << " (TrackID " << daughter.track_id
		       << " Origin " << daughter.origin
		       << " PDG " << daughter.pdg
		       << ") due to FilterNode()" << std::endl;
	  continue;
	}
	larcv::Particle part;
	try {
	  if(LArDataLabel(supera::LArDataType_t::kLArSimCh_t).empty())
	    part = MakeParticle(daughter,empty_sch_v);	  
	  else
	    part = MakeParticle(daughter,LArData<supera::LArSimCh_t>());
	}catch(const larcv::larbys& err) {
	  LARCV_NORMAL() << "Skipping a secondary (PDG,TrackID) = (" 
			 << daughter.pdg << "," << daughter.track_id << ") as it could not be turned into Particle" << std::endl;
	  continue;
	}
	if( (_filter_min_rows>0 || _filter_min_cols>0) && part.BB().empty() ) {
	  LARCV_INFO() << "Skipping a daughter " << daughter_idx
		       << " (TrackID " << daughter.track_id
		       << " Origin " << daughter.origin
		       << " PDG " << part.PdgCode() << ") due to the size of Particle" << std::endl;
	  continue;
	}
	bool skip=false;
	for(auto const& bb : part.BB()) {
	  if(bb.rows() <= _filter_min_rows && bb.cols() <= _filter_min_cols) {
	    LARCV_INFO() << "Skipping a daughter " << daughter_idx
			 << " (TrackID " << daughter.track_id
			 << " Origin " << daughter.origin
			 << " PDG " << part.PdgCode() << ") due to the size of Particle "
			 << "(row,col) = (" << bb.rows() << "," << bb.cols() << ")"
			 << std::endl;
	    skip=true;
	    break;
	  }
	}
	if(skip) continue;
	LARCV_INFO() << "Registering a daughter " << daughter_idx
		     << " (TrackID " << daughter.track_id
		     << " Origin " << daughter.origin
		     << " PDG " << part.PdgCode() << ")" << std::endl;
	sec_part_v.push_back(part);
	sec_idx_v.push_back(daughter_idx);
      }

      LARCV_INFO() << "Updating primary Particle with " << sec_part_v.size() << " children" << std::endl;
      UpdatePrimaryParticle(primary.part, sec_part_v);

      if( (_filter_min_rows>0 || _filter_min_cols>0) && primary.part.BB().empty() )
	continue;
      bool skip=false;
      for(auto const& bb : primary.part.BB()) {
	if(bb.rows() >= _filter_min_rows && bb.cols() >= _filter_min_cols)
	  continue;
	skip=true;
	break;
      }
      if(skip) continue;

      std::vector<std::pair<size_t,size_t> > part2mcnode_v;
      part2mcnode_v.clear();

      //
      // Register primary Particle
      //
      LARCV_INFO() << "Storing primary Particle (PDG " << primary.part.PdgCode()
		   << " Shape " << primary.part.Shape()
		   << " MCTIndex " << primary.part.MCTIndex()
		   << " MCSTIndex) " << primary.part.MCSTIndex()
		   << std::endl;
	
      _part_v.push_back(primary.part);

      // Record incorporated mcnode
      part2mcnode_v.emplace_back(primary_idx, larcv::kINVALID_SIZE);
      for(size_t daughter_idx=0; daughter_idx<sec_part_v.size(); ++daughter_idx) {
	auto const& secondary_idx = sec_idx_v[daughter_idx];
	auto const& secondary_part = sec_part_v[daughter_idx];
	if( (_store_g4_primary_part && secondary_part.TrackID() == secondary_part.ParentTrackID())
	    ||
	    (_store_g4_secondary_part && secondary_part.TrackID() != secondary_part.ParentTrackID())
	    )
	  part2mcnode_v.emplace_back(primary_idx,secondary_idx);
      }
      for(size_t daughter_idx=0; daughter_idx<sec_part_v.size(); ++daughter_idx) {
	auto const& daughter_part = sec_part_v[daughter_idx];
	LARCV_INFO() << "    Associated secondary (PDG " << daughter_part.PdgCode()
		     << " Shape " << daughter_part.Shape()
		     << " MCTIndex " << daughter_part.MCTIndex()
		     << " MCSTIndex) " << daughter_part.MCSTIndex()
		     << std::endl;	
      }
      _part2mcnode_vv.emplace_back(std::move(part2mcnode_v));
      
      //
      // Register secondary Particle
      //
      part2mcnode_v.clear();
      for(size_t daughter_idx=0; daughter_idx<sec_part_v.size(); ++daughter_idx) {
	auto const& secondary_idx = sec_idx_v[daughter_idx];
	auto& secondary_part = sec_part_v[daughter_idx];
	if( (_store_g4_primary_part && secondary_part.TrackID() == secondary_part.ParentTrackID())
	    ||
	    (_store_g4_secondary_part && secondary_part.TrackID() != secondary_part.ParentTrackID())
	    )
	  {
	    LARCV_INFO() << "Storing secondary Particle (PDG " << secondary_part.PdgCode()
			 << " Shape " << secondary_part.Shape()
			 << " MCTIndex " << secondary_part.MCTIndex()
			 << " MCSTIndex) " << secondary_part.MCSTIndex()
			 << std::endl;	
	    part2mcnode_v.clear();
	    part2mcnode_v.emplace_back(primary_idx,secondary_idx);
	    _part_v.emplace_back(std::move(secondary_part));
	    _part2mcnode_vv.emplace_back(std::move(part2mcnode_v));
	  }
      }
    }

    if(_store_part)
      ev_part->Emplace(std::move(_part_v));
    
    return true;
  }

  bool SuperaMCParticle::FilterNode(const supera::MCNode& node) const
  {
    if(node.source_type == supera::MCNode::SourceType_t::kMCTruth) {
      if(_pass_origin && node.origin != _pass_origin)
	return false;
    }else if(node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      auto const& mctrack = LArData<supera::LArMCTrack_t>()[node.source_index];
      LARCV_DEBUG() << "MCTrack InitE " << mctrack.Start().E() << " ... DepE "
		    << (mctrack.size() >1 ? mctrack.front().E() - mctrack.back().E() : 0) << std::endl;
      if(mctrack.Start().E() < _track_min_einit) return false;
      if(_track_min_edep > 0) {
	if(mctrack.size()<2) return false;
	if( (mctrack.front().E() - mctrack.back().E()) < _track_min_edep ) return false;
      }
      for(size_t filter_idx=0; filter_idx<_filter_pdg.size(); ++filter_idx) {
	auto const& pdg = _filter_pdg[filter_idx];
	if(pdg != mctrack.PdgCode()) continue;
	auto const& filter_min_einit = _filter_min_einit[filter_idx];
	if(mctrack.Start().E() < filter_min_einit) return false;
	auto const& filter_min_edep = _filter_min_edep[filter_idx];
	if(filter_min_edep > 0) {
	  if(mctrack.size()<2) return false;
	  if( (mctrack.front().E() - mctrack.back().E()) < filter_min_edep ) return false;
	}
      }
    }else if(node.source_type == supera::MCNode::SourceType_t::kMCShower) {
      auto const& mcshower = LArData<supera::LArMCShower_t>()[node.source_index];
      LARCV_DEBUG() << "MCShower InitE " << mcshower.Start().E() << " ... DepE "
		    << mcshower.DetProfile().E() << std::endl;
      if(mcshower.Start().E() < _shower_min_einit) return false;
      if(mcshower.DetProfile().E() < _shower_min_edep) return false;
      for(size_t filter_idx=0; filter_idx<_filter_pdg.size(); ++filter_idx) {
	auto const& pdg = _filter_pdg[filter_idx];
	if(pdg != mcshower.PdgCode()) continue;
	auto const& filter_min_einit = _filter_min_einit[filter_idx];
	if(mcshower.Start().E() < filter_min_einit) return false;
	auto const& filter_min_edep = _filter_min_edep[filter_idx];
	if( mcshower.DetProfile().E() < filter_min_edep) return false;
      }
    }
    return true;
  }

  larcv::Particle SuperaMCParticle::MakeParticle(const supera::MCNode& node,
						 const std::vector<supera::LArSimCh_t>& sch_v) const
  {
    larcv::Particle res;
    if(node.source_type == supera::MCNode::SourceType_t::kMCTruth)
      throw larbys("MCTruth type cannot make Particle using MCParticleMaker!");
    else if(node.source_type == supera::MCNode::SourceType_t::kMCTrack) {
      auto const& mctrack = LArData<supera::LArMCTrack_t>().at(node.source_index);
      if(sch_v.empty()) res = _mcpart_maker.Particle(mctrack,TimeOffset());
      else res = _mcpart_maker.Particle(mctrack,sch_v,TimeOffset());
      res.MCSTIndex(node.source_index);
    }
    else if(node.source_type == supera::MCNode::SourceType_t::kMCShower) {
      auto const& mcshower = LArData<supera::LArMCShower_t>().at(node.source_index);
      if(sch_v.empty()) res = _mcpart_maker.Particle(mcshower,TimeOffset());
      else res = _mcpart_maker.Particle(mcshower,sch_v,TimeOffset());
      res.MCSTIndex(node.source_index);
    }else
      throw larbys("Unexpected SourceType_t!");

    // format Particle
    std::vector<larcv::ImageMeta> bb_v;
    for(size_t plane=0; plane<res.BB().size(); ++plane) {
      auto const& part_meta  = res.BB().at(plane);
      auto const& event_meta = Meta().at(plane);
      bb_v.push_back(FormatMeta(part_meta,event_meta));
    }
    res.SetBB(bb_v);
    return res;
  }
  
  void SuperaMCParticle::UpdatePrimaryParticle(larcv::Particle& pri_part,
				      std::vector<larcv::Particle>& sec_part_v) const
  {
    LARCV_DEBUG() << "start" << std::endl;

    std::map<larcv::PlaneID_t, larcv::ImageMeta> sum_part_m;
    double energy_deposit = 0;
    // register primary BB
    for (auto const& bb : pri_part.BB()) {
      if(pri_part.EnergyDeposit()>0) energy_deposit += pri_part.EnergyDeposit();
      if (!(bb.height() > 1 && bb.width() > 1)) continue;
      auto iter = sum_part_m.find(bb.plane());
      if (iter == sum_part_m.end())
	sum_part_m[bb.plane()] = bb;
      else
	(*iter).second = (*iter).second.inclusive(bb);
    }

    // Next, secondary BB
    for (auto& sec_part : sec_part_v) {
      energy_deposit += sec_part.EnergyDeposit();
      // loop over plane-by-plane ImageMeta
      for (auto const& bb : sec_part.BB()) {
	if (!(bb.height() > 1 && bb.width() > 1)) continue;
	auto iter = sum_part_m.find(bb.plane());
	if (iter == sum_part_m.end())
	  sum_part_m[bb.plane()] = bb;
	else
	  (*iter).second = (*iter).second.inclusive(bb);
      }
      sec_part.MCTIndex(pri_part.MCTIndex());
    }    
    
    std::vector<larcv::ImageMeta> bb_v;
    bb_v.reserve(sum_part_m.size());

    for (auto const& plane_part : sum_part_m) {
      if(bb_v.size() <= plane_part.first)
	bb_v.resize(plane_part.first+1);
      bb_v[plane_part.first] = plane_part.second;
      LARCV_INFO() << "Updated primary Particle plane " << plane_part.first << " ... " << plane_part.second.dump();
    }
    
    pri_part.EnergyDeposit(energy_deposit);
    
    pri_part.SetBB(bb_v);
  }

  larcv::ImageMeta SuperaMCParticle::FormatMeta(const larcv::ImageMeta& part_image,
						const larcv::ImageMeta& event_image) const
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
		 << " NRows=" << rows / modular_row
		 << " NCols=" << cols / modular_col
		 << " Origin @ (" << min_x << "," << max_y << ")" << std::endl;
    larcv::ImageMeta res(width, height,
			 rows / modular_row, cols / modular_col,
			 min_x, max_y,
			 part_image.plane());

    LARCV_DEBUG() << "Event image   " << event_image.dump();

    LARCV_DEBUG() << "After format  " << res.dump();
    /*
    res = event_image.overlap(res);

    LARCV_DEBUG() << "After overlap " << res.dump();
    */
    return res;
  }
    
  void SuperaMCParticle::finalize()
  {}

}
#endif
