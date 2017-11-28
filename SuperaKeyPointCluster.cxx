#ifndef __SUPERAKEYPOINTCLUSTER_CXX__
#define __SUPERAKEYPOINTCLUSTER_CXX__

#include "SuperaKeyPointCluster.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "Voxel3DSlicer.h"

namespace larcv {

  static SuperaKeyPointClusterProcessFactory __global_SuperaKeyPointClusterProcessFactory__;

  SuperaKeyPointCluster::SuperaKeyPointCluster(const std::string name)
    : SuperaBase(name)
  {}
    
  void SuperaKeyPointCluster::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsPixel2D::configure(cfg);
    supera::ParamsVoxel3D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);

    _apply_sce    = cfg.get<bool>("ApplySCE");
    _row_pad      = cfg.get<size_t>("RowPad");
    _col_pad      = cfg.get<size_t>("ColPad");
    _pad3d_x = cfg.get<size_t>("XPad");
    _pad3d_y = cfg.get<size_t>("YPad");
    _pad3d_z = cfg.get<size_t>("ZPad");
    _cluster_primary_start   = cfg.get<bool>("UsePrimaryStart",true);
    _cluster_secondary_start = cfg.get<bool>("UseSecondaryStart",true);
    _cluster_scattering      = cfg.get<bool>("UseScattering",false);
  }

  void SuperaKeyPointCluster::initialize()
  { SuperaBase::initialize(); }

  larcv::Vertex SuperaKeyPointCluster::GetPoint(const supera::LArMCStep_t& step)
  {
    static double xyz[3];
    static const double drift_velocity = ::supera::DriftVelocity() * 1.0e-3; // make it cm/ns
    xyz[0] = step.X();
    xyz[1] = step.Y();
    xyz[2] = step.Z();
    if(_apply_sce) supera::ApplySCE(xyz);
      
    larcv::Vertex pt(xyz[0],xyz[1],xyz[2],
		     (double)(step.T() + (xyz[0] / drift_velocity)));
    
    return pt;
  }

  bool SuperaKeyPointCluster::process(IOManager& mgr)
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

    std::set<larcv::Vertex> primary_start_s;
    std::set<larcv::Vertex> secondary_start_s;
    std::set<larcv::Vertex> scattering_s;

    for(auto const& mctrack : LArData<supera::LArMCTrack_t>()) {

      if(mctrack.size()<2) continue;

      auto start = GetPoint(mctrack.Start());
      if(mctrack.TrackID() == mctrack.MotherTrackID()) {
	if(_cluster_primary_start) {
	  LARCV_INFO() << "Primary MCTrack (PDG " << mctrack.PdgCode() << ") start @ ("
		       << mctrack.Start().X() << "," << mctrack.Start().Y() << "," << mctrack.Start().Z() << ") ... w/ SCE ("
		       << start.x() << "," << start.y() << "," << start.z() << ")"
		       << " @ G4 Time = " << start.t() << " [ns]" << std::endl;
	  primary_start_s.emplace(std::move(start));
	}
      }
      else {
	if(_cluster_secondary_start) {
	  LARCV_INFO() << "Secondary MCTrack (PDG " << mctrack.PdgCode() << ") start @ ("
		       << start.x() << "," << start.y() << "," << start.z() << ")"
		       << " @ G4 Time = " << start.t() << " [ns]" << std::endl;
	  secondary_start_s.emplace(std::move(start));
	}
      }
      if(_cluster_scattering) {
	for(auto const& step : mctrack) {
	  if(step.T() == mctrack.Start().T()) continue;
	  auto pt = GetPoint(step);
	  LARCV_INFO() << "Scattering @ ("
		       << pt.x() << "," << pt.y() << "," << pt.z() << ")"
		       << " @ G4 Time = " << pt.t() << " [ns]" << std::endl;
	  scattering_s.emplace(std::move(pt));
	}
      }
    }

    for(auto const& mcshower : LArData<supera::LArMCShower_t>()) {

      if(mcshower.DetProfile().E()<=0) continue;

      auto start = GetPoint(mcshower.Start());
      
      if(mcshower.TrackID() == mcshower.MotherTrackID()) {
	if(_cluster_primary_start) {
	  LARCV_INFO() << "Primary MCShower (PDG " << mcshower.PdgCode() << ") start @ ("
		       << start.x() << "," << start.y() << "," << start.z() << ")"
		       << " @ G4 Time = " << start.t() << " [ns]" << std::endl;
	  primary_start_s.emplace(std::move(start));
	}
      }
      else {
	if(_cluster_secondary_start) {
	  LARCV_INFO() << "Primary MCShower (PDG " << mcshower.PdgCode() << ") start @ ("
		       << start.x() << "," << start.y() << "," << start.z() << ")"
		       << " @ G4 Time = " << start.t() << " [ns]" << std::endl;
	  secondary_start_s.emplace(std::move(start));
	}
      }
      if(_cluster_secondary_start) {
	if(mcshower.Start().T() != mcshower.DetProfile().T()) {
	  auto pt = GetPoint(mcshower.DetProfile());
	  if(pt == start) continue;
	  LARCV_INFO() << "DetProfile @ ("
		       << pt.x() << "," << pt.y() << "," << pt.z() << ")"
		       << " @ G4 Time = " << pt.t() << " [ns]" << std::endl;
	  secondary_start_s.emplace(std::move(pt));
	}
      }
    }

    auto meta_v = Meta();
    for(size_t plane=0; plane<meta_v.size(); ++plane) {
      
      auto& meta = meta_v[plane];
      meta.update(meta.rows() / RowCompressionFactor().at(plane),
		  meta.cols() / ColCompressionFactor().at(plane));

    }
    
    const unsigned short primary_start_val   = 3;
    const unsigned short secondary_start_val = 2;
    const unsigned short scattering_val      = 1;

    if(!(OutPixel2DLabel().empty())) {
      auto& ev_pixel2d = mgr.get_data<larcv::EventClusterPixel2D>(OutPixel2DLabel());

      for(auto const& meta : meta_v) {
	larcv::ClusterPixel2D res;
	res.meta(meta);

	LARCV_INFO() << "Creating cluster for primary start from " << primary_start_s.size() << " 3D points on ProjectionID_t " << meta.id() << std::endl;
	CreateCluster2D(res, primary_start_s, primary_start_val, TimeOffset() );

	LARCV_INFO() << "Creating cluster for secondary start from " << secondary_start_s.size() << " 3D points on ProjectionID_t " << meta.id() << std::endl;
	CreateCluster2D(res, secondary_start_s, secondary_start_val, TimeOffset() );

	LARCV_INFO() << "Creating cluster for scattering points from " << scattering_s.size() << " 3D points on ProjectionID_t " << meta.id() << std::endl;
	CreateCluster2D(res, scattering_s, scattering_val, TimeOffset() );

	ev_pixel2d.emplace(std::move(res));
      }
    }

    auto meta3d = Meta3D();
    if(!(OutVoxel3DLabel().empty())) {
      auto& ev_voxel3d = mgr.get_data<larcv::EventClusterVoxel3D>(OutVoxel3DLabel());
      ev_voxel3d.meta(meta3d);

      LARCV_INFO() << "Creating cluster for primary start from " << primary_start_s.size() << " 3D points" << std::endl;
      CreateCluster3D(ev_voxel3d, primary_start_s, primary_start_val, TimeOffset() );
      
      LARCV_INFO() << "Creating cluster for secondary start from " << secondary_start_s.size() << " 3D points" << std::endl;
      CreateCluster3D(ev_voxel3d, secondary_start_s, secondary_start_val, TimeOffset() );

      LARCV_INFO() << "Creating cluster for scattering points from " << scattering_s.size() << " 3D points" << std::endl;
      CreateCluster3D(ev_voxel3d, scattering_s, scattering_val, TimeOffset() );
    }
    
    return true;
  }

  void
  SuperaKeyPointCluster::CreateCluster3D(larcv::ClusterVoxel3D& res,
					 const std::set<larcv::Vertex>& pt_s,
					 const unsigned short val,
					 const int time_offset)
  {
    auto const& meta = res.meta();
    LARCV_INFO() << "Using Voxel3DMeta " << std::endl << meta.dump();
    std::vector<larcv::Point3D> pt_v;
    for(auto const& pt : pt_s) {
      double x = (double)(::supera::TPCG4Time2Tick(pt.t()) + time_offset + 0.5) * supera::TPCTickPeriod();
      x += supera::TriggerOffsetTPC();
      x *= ::supera::DriftVelocity();
      LARCV_INFO() << "Inspecting point (" << x << "," << pt.y() << "," << pt.z() << ")" 
		   << " ... Voxel ID: " << meta.id(x,pt.y(),pt.z()) << std::endl;
      if(meta.id(x, pt.y(), pt.z()) == larcv::kINVALID_VOXELID) 
	continue;
      larcv::Point3D xyz;
      xyz.x = x; xyz.y = pt.y(); xyz.z = pt.z();
      pt_v.emplace_back(std::move(xyz));
    }

    size_t cluster_id = res.size();
    res.resize(res.size() + pt_v.size());
    
    larcv::Point3D test_pt;

    for(size_t pt_idx=0; pt_idx<pt_v.size(); ++pt_idx) {
      auto const& pt = pt_v[pt_idx];
      auto& cluster = res.writeable_voxel_set(cluster_id);
      ++cluster_id;    

      for(size_t xpad=0; xpad<(_pad3d_x*2+1); ++xpad) {
	test_pt.x = pt.x - (xpad - _pad3d_x) * meta.size_voxel_x();
	if(test_pt.x < meta.min_x() || meta.max_x() < test_pt.x) continue;

	for(size_t ypad=0; ypad<(_pad3d_y*2+1); ++ypad) {
	  test_pt.y = pt.y - (ypad - _pad3d_y) * meta.size_voxel_y();
	  if(test_pt.y < meta.min_y() || meta.max_y() < test_pt.y) continue;

	  for(size_t zpad=0; zpad<(_pad3d_z*2+1); ++zpad) {
	    test_pt.z = pt.z - (zpad - _pad3d_z) * meta.size_voxel_z();
	    if(test_pt.z < meta.min_z() || meta.max_z() < test_pt.z) continue;

	    auto voxel_id = meta.id(test_pt);
	    if(voxel_id == larcv::kINVALID_VOXELID) continue;
	    
	    larcv::Voxel vox(voxel_id, (float)val);
	    
	    cluster.emplace(std::move(vox),true);

	  }
	}
      }
    }
    LARCV_INFO() << "Clustered " << res.size() << " points" << std::endl;
  }

  void
  SuperaKeyPointCluster::CreateCluster2D(larcv::ClusterPixel2D& res,
					 const std::set<larcv::Vertex>& pt_s,
					 const unsigned short val,
					 const int time_offset)
  {

    double xyz[3];
    auto const& meta = res.meta();
    LARCV_DEBUG() << "Clustering ProjectionID_t " << meta.id() << std::endl;
    std::vector<std::pair<int,int> > wire_tick_v;

    for(auto const& pt : pt_s) {
      xyz[0] = pt.x();
      xyz[1] = pt.y();
      xyz[2] = pt.z();
      double tick_float = (double)(::supera::TPCG4Time2Tick(pt.t()) + time_offset) + supera::PlaneTickOffset(0,meta.id()) + 0.5;
      double wire_float = (double)(::supera::NearestWire(xyz, meta.id())) + 0.5;
      LARCV_DEBUG() << wire_float << " ... " << meta.min_x() << " => " << meta.max_x() << std::endl;
      LARCV_DEBUG() << tick_float << " ... " << meta.min_y() << " => " << meta.max_y() << std::endl;
      if(wire_float < meta.min_x() || meta.max_x() < wire_float) continue;
      if(tick_float < meta.min_y() || meta.max_y() < tick_float) continue;
      
      auto tick = meta.row(tick_float);
      auto wire = meta.col(wire_float);
      wire_tick_v.emplace_back((int)wire,(int)tick);
    }

    size_t cluster_id = res.size();
    res.resize(res.size()+wire_tick_v.size());
    for(size_t pt_idx=0; pt_idx<wire_tick_v.size(); ++pt_idx) {

      auto const& wire = wire_tick_v[pt_idx].first;
      auto const& tick = wire_tick_v[pt_idx].second;
      auto& cluster = res.writeable_voxel_set(cluster_id);
      ++cluster_id;

      LARCV_DEBUG() << "Masking around (tick,wire) = (" << tick << "," << wire << ")" << std::endl;
	
      for(int row=((int)tick - (int)_row_pad); row<=((int)tick + (int)(_row_pad)); ++row) {
	for(int col=((int)wire - (int)_col_pad); col<=((int)wire + (int)_col_pad); ++col) {
	  if(row<0 || row >= (int)(meta.rows())) continue;
	  if(col<0 || col >= (int)(meta.cols())) continue;
	  LARCV_DEBUG() << "    Registering (row,col) = (" << row << "," << col << ") w/ value " << val << std::endl;
	  auto voxel_id = meta.index((size_t)row,(size_t)col);
	  larcv::Voxel vox(voxel_id,(float)val);
	  cluster.emplace(std::move(vox),true);
	}
      }
    }
    LARCV_INFO() << "Clustered " << res.size() << " points @ projection " << meta.id() << std::endl;

  }
  
  void SuperaKeyPointCluster::finalize()
  { SuperaBase::finalize(); }

}
#endif
