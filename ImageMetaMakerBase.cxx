#ifndef __IMAGEMETAMAKERBASE_CXX__
#define __IMAGEMETAMAKERBASE_CXX__

#include "ImageMetaMakerBase.h"

namespace supera {

  void ImageMetaMakerBase::configure(const Config_t& cfg)
  {
    _comp_rows  = cfg.get<std::vector<size_t> >("EventCompRows");
    _comp_cols  = cfg.get<std::vector<size_t> >("EventCompCols");

    if(_comp_rows.size() != _comp_cols.size()) {
      std::cerr << "EventCompRows size != EventCompCols size" << std::endl;
      throw std::exception();
    }
  }  
}
#endif
