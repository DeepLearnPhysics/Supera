#ifndef __SUPERA_LAR2IMAGE_H__
#define __SUPERA_LAR2IMAGE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "larcv/core/DataFormat/DataFormatTypes.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/Voxel2D.h"
#include "larcv/core/DataFormat/Voxel3D.h"

namespace supera {

  //
  // Hit => Image2D
  //
  larcv::Image2D Hit2Image2D(const larcv::ImageMeta& meta,
			     const std::vector<supera::LArHit_t>& hits,
			     const int time_offset=0,
			     const bool smear=false);
  
  std::vector<larcv::Image2D> Hit2Image2D(const std::vector<larcv::ImageMeta>& meta_v,
					  const std::vector<supera::LArHit_t>& hits,
					  const int time_offset=0,
					  const bool smear=false);
  
  //
  // Wire => Image2D
  //
  larcv::Image2D Wire2Image2D(const larcv::ImageMeta& meta,
			      const std::vector<supera::LArWire_t>& wires,
			      const int time_offset=0);
  
  std::vector<larcv::Image2D> Wire2Image2D(const std::vector<larcv::ImageMeta>& meta_v,
					   const std::vector<supera::LArWire_t>& wires,
					   const int time_offset=0);

  //
  // OpDigit => Image2D
  //
  larcv::Image2D OpDigit2Image2D(const larcv::ImageMeta& meta,
				 const std::vector<supera::LArOpDigit_t>& opdigit_v,
				 int time_offset=0);

  //
  // SimChannel => Image2D
  //
  /*
  std::vector<larcv::Image2D> SimCh2Image2D(const std::vector<larcv::ImageMeta>& meta_v,
					    const std::vector<larcv::ROIType_t>& track2type_v,
					    const std::vector<supera::LArSimCh_t>& sch_v,
					    const int time_offset);
  */
  //
  // SimChannel => Voxel3D
  //
  void
  SimCh2ClusterVoxel3D(larcv::ClusterVoxel3D& res,
                       const std::vector<supera::LArSimCh_t>& sch_v,
                       const std::vector<size_t>& trackid2cluster,
                       const int time_offset,
		       const bool use_true_pos,
                       const larcv::ProjectionID_t id = larcv::kINVALID_PROJECTIONID);

  void
  SimCh2SparseTensor3D(larcv::SparseTensor3D& res,
		       const std::vector<supera::LArSimCh_t>& sch_v,
		       const std::vector<size_t>& track_v,
		       const int time_offset,
		       const bool use_true_pos,
		       const larcv::ProjectionID_t id = larcv::kINVALID_PROJECTIONID);

  //
  // SimChannel => Pixel2DCluster
  //
  void
  SimCh2ClusterPixel2D(std::vector<larcv::ClusterPixel2D>& res,
		       const std::vector<supera::LArSimCh_t>& sch_v,
		       const std::vector<size_t>& trackid2cluster,
		       const int time_offset);

}
#endif
//#endif
//#endif
