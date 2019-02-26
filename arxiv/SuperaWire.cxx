#ifndef __SUPERAWIRE_CXX__
#define __SUPERAWIRE_CXX__

#include "SuperaWire.h"
#include "LAr2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "Voxel3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"

namespace larcv {

  static SuperaWireProcessFactory __global_SuperaWireProcessFactory__;

  SuperaWire::SuperaWire(const std::string name)
    : SuperaBase(name)
  {}

  void SuperaWire::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
  }

  void SuperaWire::initialize()
  { SuperaBase::initialize(); }

  bool SuperaWire::process(IOManager& mgr)
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

    auto image_v = supera::Wire2Image2D(meta_v, LArData<supera::LArWire_t>(), TimeOffset());
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

  void SuperaWire::finalize()
  {}


}
#endif
