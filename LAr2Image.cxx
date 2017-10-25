#ifndef __SUPERA_LAR2IMAGE_CXX__
#define __SUPERA_LAR2IMAGE_CXX__

#include "LAr2Image.h"
#include "larcv/core/Base/larcv_logger.h"

namespace supera {

  larcv::Image2D Hit2Image2D(const larcv::ImageMeta & meta,
                             const std::vector<supera::LArHit_t>& hits, const int time_offset)
  {
    //int nticks = meta.rows();
    //int nwires = meta.cols();
    //size_t row_comp_factor = (size_t)(meta.pixel_height());
    const int ymax = meta.max_y() - 1; // Need in terms of row coordinate
    const int ymin = (meta.min_y() >= 0 ? meta.min_y() : 0);
    larcv::Image2D img(meta);

    LARCV_SINFO() << "Filling an image: " << meta.dump();
    LARCV_SINFO() << "(ymin,ymax) = (" << ymin << "," << ymax << ")" << std::endl;

    for (auto const& h : hits) {
      auto const projection_id = ::supera::ChannelToProjectionID(h.Channel());
      auto const image_x = ::supera::ChannelToImageX(h.Channel());

      if (projection_id != meta.id()) continue;

      size_t col = 0;
      try { col = meta.col(image_x); }
      catch (const larcv::larbys&) { continue; }

      int row = int(h.PeakTime() + 0.5) + time_offset;
      if (row > ymax || row < ymin) continue;
      img.set_pixel(ymax - row, col, h.Integral());
    }
    return img;
  }

  std::vector<larcv::Image2D> Hit2Image2D(const std::vector<larcv::ImageMeta>& meta_v,
                                          const std::vector<supera::LArHit_t>& hits,
                                          const int time_offset)
  {
    std::vector<larcv::Image2D> res_v;
    for (size_t p = 0; p < ::supera::Nplanes(); ++p) {
      auto const& meta = meta_v.at(p);
      res_v.emplace_back(std::move(Hit2Image2D(meta, hits, time_offset)));
      res_v.back().index(res_v.size() - 1);
    }
    return res_v;
  }

  larcv::Image2D Wire2Image2D(const larcv::ImageMeta& meta,
                              const std::vector<supera::LArWire_t>& wires,
                              const int time_offset)
  {
    //int nticks = meta.rows();
    //int nwires = meta.cols();
    size_t row_comp_factor = (size_t)(meta.pixel_height());
    const int ymax = meta.max_y() - 1; // Need in terms of row coordinate
    const int ymin = (meta.min_y() >= 0 ? meta.min_y() : 0);
    larcv::Image2D img(meta);
    img.paint(0.);

    LARCV_SINFO() << "Filling an image: " << meta.dump();
    LARCV_SINFO() << "(ymin,ymax) = (" << ymin << "," << ymax << ")" << std::endl;

    for (auto const& wire : wires) {
      auto const projection_id = ::supera::ChannelToProjectionID(wire.Channel());
      auto const image_x = ::supera::ChannelToImageX(wire.Channel());

      if (projection_id != meta.id()) continue;

      size_t col = 0;
      try {
        col = meta.col(image_x);
      } catch (const ::larcv::larbys&) {
        continue;
      }

      for (auto const& range : wire.SignalROI().get_ranges()) {

        auto const& adcs = range.data();
        //double sumq = 0;
        //for(auto const& v : adcs) sumq += v;
        //sumq /= (double)(adcs.size());
        //if(sumq<3) continue;

        int start_index = range.begin_index() + time_offset;
        int end_index   = start_index + adcs.size() - 1;
        if (start_index > ymax || end_index < ymin) continue;

        if (row_comp_factor > 1) {

          for (size_t index = 0; index < adcs.size(); ++index) {
            if ((int)index + start_index < ymin) continue;
            if ((int)index + start_index > ymax) break;
            auto row = meta.row((double)(start_index + index));
            img.set_pixel(row, col, adcs[index] + img.pixel(row, col));
          }
        } else {
          // Fill matrix from start_index => end_index of matrix row
          // By default use index 0=>length-1 index of source vector
          int nskip = 0;
          int nsample = adcs.size();
          if (end_index   > ymax) {
            LARCV_SDEBUG() << "End index (" << end_index << ") exceeding image bound (" << ymax << ")" << std::endl;
            nsample   = adcs.size() - (end_index - ymax);
            end_index = ymax;
            LARCV_SDEBUG() << "Corrected End index = " << end_index << std::endl;
          }
          if (start_index < ymin) {
            LARCV_SDEBUG() << "Start index (" << start_index << ") exceeding image bound (" << ymin << ")" << std::endl;
            nskip = ymin - start_index;
            nsample -= nskip;
            start_index = ymin;
            LARCV_SDEBUG() << "Corrected Start index = " << start_index << std::endl;
          }
          LARCV_SDEBUG() << "Calling Image2D::copy..." << std::endl
                         << "      source wf : start index = " << range.begin_index() << " length = " << adcs.size() << std::endl
                         << "      (row,col) : (" << start_index << "," << col << ")" << std::endl
                         << "      nskip     : "  << nskip << std::endl
                         << "      nsample   : "  << nsample << std::endl;
          try {
            //img.reverse_copy(ymax - end_index,
            img.copy(meta.row(start_index),
                     col,
                     adcs,
                     nsample);
          } catch (const ::larcv::larbys& err) {
            LARCV_SCRITICAL() << "Attempted to fill an image..." << std::endl
                              << meta.dump()
                              << "(ymin,ymax) = (" << ymin << "," << ymax << ")" << std::endl
                              << "Called a reverse_copy..." << std::endl
                              << "      source wf : plane = " << projection_id << " wire = " << image_x << std::endl
                              << "      timing    : start index = " << range.begin_index() << " length = " << adcs.size() << std::endl
                              << "      (row,col) : (" << start_index << "," << col << ")" << std::endl
                              << "      nskip     : "  << nskip << std::endl
                              << "Re-throwing an error:" << std::endl;
            throw err;
          }
        }
      }
    }
    return img;
  }

  std::vector<larcv::Image2D>
  Wire2Image2D(const std::vector<larcv::ImageMeta>& meta_v,
               const std::vector<supera::LArWire_t>& wires,
               const int time_offset)
  {
    std::vector<larcv::Image2D> res_v;
    for (size_t p = 0; p < ::supera::Nplanes(); ++p) {
      auto const& meta = meta_v.at(p);
      res_v.emplace_back(std::move(Wire2Image2D(meta, wires, time_offset)));
      res_v.back().index(res_v.size() - 1);
    }
    return res_v;
  }

  larcv::Image2D OpDigit2Image2D(const larcv::ImageMeta& meta,
                                 const std::vector<supera::LArOpDigit_t>& opdigit_v,
                                 int time_offset)
  {
    larcv::Image2D img(meta);
    img.paint(2048);

    std::vector<float> tmp_wf(meta.rows(), 2048);
    for (auto const& opdigit : opdigit_v) {
      if (opdigit.size() < 1000) continue;
      auto const col = opdigit.ChannelNumber();
      if (meta.min_x() > col) continue;
      if (col >= meta.max_x()) continue;
      //
      // HACK: right way is to use TimeService + trigger.
      //       for now I just record PMT beamgate tick=0 as edge of an image (w/ offset)
      //
      size_t nskip = 0;
      if (time_offset < 0) nskip = (-1 * time_offset);
      if (nskip >= opdigit.size()) continue;
      for (auto& v : tmp_wf) v = 2048;
      size_t num_pixel = std::min(meta.rows(), opdigit.size() - nskip);
      for (size_t i = 0; i < num_pixel; ++i) tmp_wf[i] = (float)(opdigit[nskip + i]);
      img.copy(0, col, &(tmp_wf[0]), num_pixel);
      //img.reverse_copy(0,col,opdigit,nskip,num_pixel);
    }
    return img;
  }
  /*
  std::vector<larcv::Image2D>
  SimCh2Image2D(const std::vector<larcv::ImageMeta>& meta_v,
    const std::vector<larcv::ROIType_t>& track2type_v,
    const std::vector<supera::LArSimCh_t>& sch_v,
    const int time_offset)
  {
    LARCV_SINFO() << "Filling semantic-segmentation ground truth image..." << std::endl;
    std::vector<larcv::Image2D> img_v;
    for (auto const& meta : meta_v) {
      LARCV_SINFO() << meta.dump();
      img_v.emplace_back(larcv::Image2D(meta));
    }

    static std::vector<float> column;
    for (auto const& img : img_v) {
      if (img.meta().rows() >= column.size())
  column.resize(img.meta().rows() + 1, (float)(::larcv::kROIUnknown));
    }

    for (auto const& sch : sch_v) {
      auto ch = sch.Channel();
      auto const& wid = ::supera::ChannelToWireID(ch);
      auto const& plane = wid.Plane;
      auto& img = img_v.at(plane);
      auto const& meta = img.meta();

      size_t col = wid.Wire;
      if (col < meta.min_x()) continue;
      if (meta.max_x() <= col) continue;
      if (plane != img.meta().id()) continue;

      col -= (size_t)(meta.min_x());

      // Initialize column vector
      for (auto& v : column) v = (float)(::larcv::kROIUnknown);
      //for (auto& v : column) v = (float)(-1);

      for (auto const tick_ides : sch.TDCIDEMap()) {
  int tick = supera::TPCTDC2Tick((double)(tick_ides.first)) + time_offset;
  if (tick < meta.min_y()) continue;
  if (tick >= meta.max_y()) continue;
  // Where is this tick in column vector?
  size_t index = (size_t)(meta.max_y() - tick);
  // Pick type
  double energy = 0;
  ::larcv::ROIType_t roi_type =::larcv::kROIUnknown;
  for (auto const& edep : tick_ides.second) {
    if (edep.energy < energy) continue;
    if (std::abs(edep.trackID) >= (int)(track2type_v.size())) continue;
    auto temp_roi_type = track2type_v[std::abs(edep.trackID)];
    if (temp_roi_type ==::larcv::kROIUnknown) continue;
    energy = edep.energy;
    roi_type = (::larcv::ROIType_t)temp_roi_type;
  }
  column[index] = roi_type;
      }
      // mem-copy column vector
      img.copy(0, col, column, img.meta().rows());
    }
    return img_v;
  }
  */

  larcv::Voxel3DSet
  SimCh2Voxel3D(const larcv::Voxel3DMeta& meta,
    const std::vector<int>& track_v,
    const std::vector<supera::LArSimCh_t>& sch_v,
    const int time_offset,
    const size_t plane)
  {
    LARCV_SINFO() << "Filling Voxel3D ground truth volume..." << std::endl;
    larcv::Voxel3DSet res(meta);
    //double x, y, z, x_tick;
    double y, z, x_tick;
    //std::cout << "x_offset " << x_offset << std::endl;
    for (auto const& sch : sch_v) {
      auto ch = sch.Channel();

      auto const& wid = ::supera::ChannelToWireID(ch);
      if( plane != (size_t)(wid.Plane) ) continue;

      for (auto const tick_ides : sch.TDCIDEMap()) {

  x_tick = (supera::TPCTDC2Tick(tick_ides.first) + time_offset) * supera::TPCTickPeriod()  * supera::DriftVelocity();

  //std::cout << tick_ides.first << " : " << supera::TPCTDC2Tick(tick_ides.first) << " : " << time_offset << " : " << x_tick << std::flush;
  //std::cout << (supera::TPCTDC2Tick(tick_ides.first) + time_offset) << std::endl
  //<< (supera::TPCTDC2Tick(tick_ides.first) + time_offset) * supera::TPCTickPeriod() << std::endl
  //<< (supera::TPCTDC2Tick(tick_ides.first) + time_offset) * supera::TPCTickPeriod() * supera::DriftVelocity() << std::endl
  //<< std::endl;

  for (auto const& edep : tick_ides.second) {

    if (std::abs(edep.trackID) >= (int)(track_v.size())) continue;

    if (track_v[std::abs(edep.trackID)]<=0) continue;

    //x = edep.x;
    y = edep.y;
    z = edep.z;
    //supera::ApplySCE(x,y,z);
    //std::cout << " ... " << x << std::endl;
    // Now use tick-based position for x
    auto vid = meta.ID(x_tick,y,z);
    if(vid == larcv::kINVALID_VOXEL3DID) continue;

    larcv::Voxel3D vx(vid,edep.energy);
    res.Emplace(std::move(vx));
  }
      }
    }
    return res;
  }

  larcv::VoxelSetArray
  SimCh2Voxel3DCluster(const larcv::Voxel3DMeta& meta,
                       const std::vector<supera::LArSimCh_t>& sch_v,
                       const std::vector<size_t>& trackid2cluster,
                       const int time_offset)
  {

    // Create the data component of VoxelSetArray2D
    // Note: voxel_vv is VoxelSet (particle) array
    larcv::VoxelSetArray voxel_vv;
    voxel_vv.meta(meta);
    // figure out # of clusters to be made
    size_t num_clusters = 0;
    for (auto const& cidx : trackid2cluster) {
      if (cidx == larcv::kINVALID_SIZE) continue;
      if (cidx >= num_clusters) num_clusters = cidx + 1;
    }
    voxel_vv.resize(num_clusters+1);

    // Loop over sim channels
    for (auto const& sch : sch_v) {
      auto ch = sch.Channel();
      // Loop over ticks
      for (auto const tick_ides : sch.TDCIDEMap()) {
	larcv::Voxel vox;
	// Loop over energy deps!
	for (auto const& edep : tick_ides.second) {
	  if(edep.energy <= 0) continue;
	  // figure out cluster id
	  size_t vox_id = meta.id(edep.x,edep.y,edep.z);
	  if(vox_id == kINVALID_VOXELID) continue;
	  size_t trackid = std::abs(edep.trackID);
	  size_t cluster_id = num_clusters;
	  if (trackid < trackid2cluster.size() && 
	      trackid2cluster[trackid] != larcv::kINVALID_SIZE)
	    cluster_id = trackid2cluster[trackid];
	  // Fill voxel
	  vox.set(vox_id,edep.energy);
	  auto& voxel_v = voxel_vv.writeable_voxel_set(cluster_id);
	  voxel_v.add(vox);
	}
      }
    }
  
    return voxel_vv;
  }

  larcv::VoxelSetArray2D
  SimCh2Pixel2DCluster(const std::vector<larcv::ImageMeta>& meta_v,
                       const std::vector<supera::LArSimCh_t>& sch_v,
                       const std::vector<size_t>& trackid2cluster,
                       const int time_offset)
  {

    // Create the data component of VoxelSetArray2D
    // Note: voxel_vvv is VoxelSet (particle) array per projection id
    std::vector<larcv::VoxelSetArray> voxel_vvv;
    // figure out # of clusters to be made
    size_t num_clusters = 0;
    for (auto const& cidx : trackid2cluster) {
      if (cidx == larcv::kINVALID_SIZE) continue;
      if (cidx >= num_clusters) num_clusters = cidx + 1;
    }
    // Create projection id => meta mapping, also set size of clusters-per-plane in voxel_vvv
    std::vector<size_t> projection_id_to_meta_index;
    for(size_t i=0; i<meta_v.size(); ++i) {
      auto const& meta = meta_v[i];
      if(projection_id_to_meta_index.size() <= meta.id()) {
        projection_id_to_meta_index.resize(meta.id()+1,larcv::kINVALID_SIZE);
	voxel_vvv.resize(meta.id()+1);
      }
      voxel_vvv[meta.id()].resize(num_clusters+1); // 1 bigger to contain "unknown"
      projection_id_to_meta_index[meta.id()] = i;
    }
    // Loop over sim channels
    for (auto const& sch : sch_v) {
      auto ch = sch.Channel();
      // Figure out image-x-coordinate from channel
      double x = ::supera::ChannelToImageX(ch);
      // Figure out meta projection id
      larcv::ProjectionID_t projection_id = supera::ChannelToProjectionID(ch);
      // Get meta if found. else continue
      if(projection_id_to_meta_index.size() <= projection_id) continue;
      if(projection_id_to_meta_index[projection_id] == larcv::kINVALID_SIZE) continue;
      auto const& meta = meta_v[projection_id_to_meta_index[projection_id]];
      // Figure out column
      if(x < meta.min_x() || meta.max_x() <= x) continue;
      size_t col = meta.col(x);
      // Get VoxelSetArray of this projection
      auto& voxel_vv = voxel_vvv[projection_id];
      // Loop over ticks
      for (auto const tick_ides : sch.TDCIDEMap()) {
	// Apply time offset, continue if out-of-range of image meta
	double tick = supera::TPCTDC2Tick((double)(tick_ides.first)) + time_offset;
	if (tick < meta.min_y() || meta.max_y() <= tick) continue;
	size_t row = meta.row(tick);
	// Create temporary voxel
	larcv::Voxel vox(meta.index(row,col), 0.);
	// Loop over energy deps!
	for (auto const& edep : tick_ides.second) {
	  if(edep.energy <= 0) continue;
	  // figure out cluster id
	  size_t trackid = std::abs(edep.trackID);
	  size_t cluster_id = num_clusters;
	  if (trackid < trackid2cluster.size() && 
	      trackid2cluster[trackid] != larcv::kINVALID_SIZE)
	    cluster_id = trackid2cluster[trackid];
	  // Fill voxel
	  vox.set(vox.id(),edep.energy);
	  auto& voxel_v = voxel_vv.writeable_voxel_set(cluster_id);
	  voxel_v.add(vox);
	}
      }
    }
  
    larcv::VoxelSetArray2D res;
    for(auto const& meta : meta_v) {
      auto meta_copy = meta;
      res.emplace(std::move(voxel_vvv[meta_copy.id()]),std::move(meta_copy));
    }
    return res;
  }


}
#endif
