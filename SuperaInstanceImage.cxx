#ifndef __SUPERALARFLOW_CXX__
#define __SUPERALARFLOW_CXX__

#include "SuperaInstanceImage.h"
#include "Instance2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"
//#include "larcv/core/DataFormat/DataFormatUtil.h"

namespace larcv {

  static SuperaInstanceImageProcessFactory __global_SuperaInstanceImageProcessFactory__;
  
  SuperaInstanceImage::SuperaInstanceImage(const std::string name)
    : SuperaBase(name)
  {}
  
  void SuperaInstanceImage::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    _origin = cfg.get<unsigned short>("Origin",0);
    m_ancestor_label = cfg.get<std::string>("AncestorImageLabel");
  }

  void SuperaInstanceImage::initialize()
  {}

  bool SuperaInstanceImage::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }
    
    // get meta for image we initially fill
    auto const& filler_meta_v = Meta();
    // get meta of output image (the shared meta)
    ImageMetaMaker tmp_maker;
    auto const& shared_meta_v = tmp_maker.Meta(); // since tmp_maker is null, will return shared meta
    
    if(filler_meta_v.empty()) {
      LARCV_CRITICAL() << "Filler meta not created!" << std::endl;
      throw larbys();
    }
    if(shared_meta_v.empty()) {
      LARCV_CRITICAL() << "Shared meta not created!" << std::endl;
      throw larbys();
    }

    auto ev_image = (EventImage2D*)(mgr.get_data("image2d",OutImageLabel()));
    if(!ev_image) {
      LARCV_CRITICAL() << "Output image could not be created!" << std::endl;
      throw larbys();
    }
    if(!(ev_image->image2d_array().empty())) {
      LARCV_CRITICAL() << "Output image array not empty!" << std::endl;
      throw larbys();
    }
    
    auto ev_ancestor = (EventImage2D*)(mgr.get_data("image2d",m_ancestor_label));
    if(!ev_ancestor) {
      LARCV_CRITICAL() << "Ancestor image could not be created!" << std::endl;
      throw larbys();
    }
    if(!(ev_ancestor->image2d_array().empty())) {
      LARCV_CRITICAL() << "Ancestor image array not empty!" << std::endl;
      throw larbys();
    }


    // the map I need to make
    // [trackid] -> [ancenstor id]

    LARCV_SDEBUG() << "==============================================" << std::endl;
    LARCV_SDEBUG() << "MC Track Scraping" << std::endl;

    std::map<int,int> trackid2ancestorid;

    for(auto const& mctrack : LArData<supera::LArMCTrack_t>()) {
      // std::cout << "mctrack: "
      // 		<< " id=" << mctrack.TrackID() 
      // 		<< " ancestorid=" << mctrack.AncestorTrackID() 
      // 		<< " motherid=" << mctrack.MotherTrackID() 
      // 		<< " pdg=" << mctrack.PdgCode() 
      // 		<< " origin=" << mctrack.Origin() 
      // 		<< std::endl;

      if(_origin && ((unsigned short)(mctrack.Origin())) != _origin) continue;
      
      trackid2ancestorid[mctrack.TrackID()] = mctrack.AncestorTrackID();
    }

    LARCV_SDEBUG() << "==============================================" << std::endl;
    LARCV_SDEBUG() << "MC Shower Scraping" << std::endl;
    for(auto const& mcshower : LArData<supera::LArMCShower_t>()) {

      // std::cout << "mcshower: "
      // 		<< " id=" << mcshower.TrackID() 
      // 		<< " ancestorid=" << mcshower.AncestorTrackID() 
      // 		<< " motherid=" << mcshower.MotherTrackID() 
      // 		<< " pdg=" << mcshower.PdgCode() 
      // 		<< " origin=" << mcshower.Origin() 
      // 		<< std::endl;

      if(_origin && ((unsigned short)(mcshower.Origin())) != _origin) continue;
      
      trackid2ancestorid[mcshower.TrackID()] = mcshower.AncestorTrackID();
    }


    // std::vector<float> row_compression_factor;
    // std::vector<float> col_compression_factor;
    // for ( auto const& meta : shared_meta_v ) {
    //   row_compression_factor.push_back( RowCompressionFactor().at(meta.id()) );
    //   col_compression_factor.push_back( ColCompressionFactor().at(meta.id()) );
    // }

    std::vector<larcv::Image2D> idimg_v;
    std::vector<larcv::Image2D> ancestor_v;
    supera::Instance2Image(filler_meta_v, shared_meta_v, 
			   trackid2ancestorid, LArData<supera::LArSimCh_t>(),
			   TimeOffset(),
			   idimg_v, ancestor_v );
          
    ev_image->emplace(std::move(idimg_v));
    ev_ancestor->emplace(std::move(ancestor_v));
    
    return true;
  }

  void SuperaInstanceImage::finalize()
  {}

}
#endif
