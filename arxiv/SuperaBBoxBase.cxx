#ifndef __SUPERAMETAMAKER_CXX__
#define __SUPERAMETAMAKER_CXX__

#include "SuperaBBoxBase.h"
#include "SuperaCSVReader.h"
#include "LAr2Image.h"
#include "ImageBBoxBase.h"
#include "PulledPork3DSlicer.h"
#include "Voxel3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"

namespace larcv {

  static SuperaBBoxBaseProcessFactory __global_SuperaBBoxBaseProcessFactory__;

  SuperaBBoxBase::SuperaBBoxBase(const std::string name)
    : SuperaBase(name)
    , _meta_maker(nullptr)
  {}
  
  bool SuperaBBoxBase::is(const std::string question) const
  {
    if(question == "Supera") return true;
    if(question == "SuperaBBoxBase") return true;
    return false;
  }

  void SuperaBBoxBase::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    if(_meta_maker) delete _meta_maker;
    _meta_maker = supera::CreateImageBBoxBase(cfg);
  }

  void SuperaBBoxBase::initialize()
  { SuperaBase::initialize(); }

  bool SuperaBBoxBase::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    return true;
  }

  void SuperaBBoxBase::finalize()
  {}


}
#endif
