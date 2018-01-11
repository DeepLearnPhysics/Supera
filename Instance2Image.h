#ifndef __SUPERA_INSTANCE_LAR2IMAGE_H__
#define __SUPERA_INSTANCE_LAR2IMAGE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include <map>

namespace supera {

  //
  // SimChannel => Instance/Ancestor Images
  // 
  void Instance2Image( const std::vector<larcv::ImageMeta>& meta_v,
		       const std::vector<larcv::ImageMeta>& out_meta_v,
		       const std::map<int,int>& trackid2ancestorid,
		       const std::vector<supera::LArSimCh_t>& sch_v,
		       const int time_offset,
		       std::vector<larcv::Image2D>& img_out_v,
		       std::vector<larcv::Image2D>& ancestor_out_v );
  
}
#endif
//#endif
//#endif
