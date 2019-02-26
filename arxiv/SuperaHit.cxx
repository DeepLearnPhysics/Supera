#ifndef __SUPERAHIT_CXX__
#define __SUPERAHIT_CXX__

#include "SuperaHit.h"
#include "LAr2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "Voxel3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"
namespace larcv {

  static SuperaHitProcessFactory __global_SuperaHitProcessFactory__;

  SuperaHit::SuperaHit(const std::string name)
    : SuperaBase(name)
  {}
    
  void SuperaHit::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    _smear = cfg.get<bool>("Smear",false);
  }

  void SuperaHit::initialize()
  { SuperaBase::initialize(); }

  bool SuperaHit::process(IOManager& mgr)
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

    auto const& meta_v = Meta();

    if(meta_v.empty()) {
      LARCV_CRITICAL() << "Meta not created!" << std::endl;
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

    auto image_v = supera::Hit2Image2D(meta_v, LArData<supera::LArHit_t>(), TimeOffset(), _smear);
    /*
    for(size_t plane=0; plane<image_v.size(); ++plane) {
      auto& image = image_v[plane];
      image.compress(image.meta().rows() / RowCompressionFactor().at(plane),
		     image.meta().cols() / ColCompressionFactor().at(plane),
		     larcv::Image2D::kSum);
    }
    */
    ev_image->emplace(std::move(image_v));
    
    return true;
  }

  void SuperaHit::finalize()
  {SuperaBase::finalize();}

}
#endif
