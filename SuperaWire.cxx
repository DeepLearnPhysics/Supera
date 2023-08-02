#ifndef __SUPERAWIRE_CXX__
#define __SUPERAWIRE_CXX__

#include "SuperaWire.h"
#include "larcv/core/DataFormat/EventVoxel2D.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaWireProcessFactory __global_SuperaWireProcessFactory__;

  SuperaWire::SuperaWire(const std::string name)
    : SuperaBase(name)
  {}

  int SuperaWire::plane_index(unsigned int cryo_id, unsigned int tpc_id, unsigned int plane_id)
  {
    if(_scan.size() <= cryo_id)
      return -1;
    if(_scan[cryo_id].size() <= tpc_id)
      return -1;
    if(_scan[cryo_id][tpc_id].size() <= plane_id)
      return -1;

    return _scan[cryo_id][tpc_id][plane_id];
  }

  void SuperaWire::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    //_image_time_ticks = cfg.get<size_t>      ("ImageTimeTicks",4);
    _npx_rows = cfg.get<int>("RowCount");
    _npx_columns = cfg.get<int>("ColumnCount");
    _time_compression = cfg.get<int>("TimeCompression");
    assert(_npx_rows>0 && _npx_columns>0 && _time_compression>0);
    _output_producer  = cfg.get<std::string> ("OutputProducer");

    _ref_meta3d_cluster3d = cfg.get<std::string>("Meta3DFromCluster3D", "");
    _ref_meta3d_tensor3d  = cfg.get<std::string>("Meta3DFromTensor3D",  "");
    
    _adc_threshold = cfg.get<double>("ADCThreshold",0);
    // construct _scan as 3d array of [#cryo][#tpc][#plane]
    auto cryostat_v   = cfg.get<std::vector<unsigned short> >("CryostatList");
    auto tpc_v        = cfg.get<std::vector<unsigned short> >("TPCList");
    auto plane_v      = cfg.get<std::vector<unsigned short> >("PlaneList");

    auto geop = lar::providerFrom<geo::Geometry>();
    _scan.resize(geop->Ncryostats());
    for(size_t cryoid=0; cryoid<_scan.size(); ++cryoid) {
      auto const& cryostat = geop->Cryostat(geo::CryostatID(cryoid));
      auto& scan_cryo = _scan[cryoid];
      scan_cryo.resize(cryostat.NTPC());
      for(size_t tpcid=0; tpcid<scan_cryo.size(); ++tpcid) {
	auto const& tpc = cryostat.TPC(tpcid);
	auto& scan_tpc = scan_cryo[tpcid];
	scan_tpc.resize(tpc.Nplanes(),-1);
      }
    }

    _valid_nplanes = 0;
    for(auto const& cryo_id : cryostat_v) {
      auto const& cryostat = geop->Cryostat(geo::CryostatID(cryo_id));
      if(_scan.size()<=cryo_id)
	_scan.resize(cryo_id+1);
      for(auto const& tpc_id : tpc_v) {
	if(!cryostat.HasTPC(tpc_id)) continue;
	auto const& tpc = cryostat.TPC(tpc_id);
	if(_scan[cryo_id].size() <= tpc_id)
	  _scan[cryo_id].resize(tpc_id+1);
	for(auto const& plane_id : plane_v) {
	  if(!tpc.HasPlane(plane_id)) continue;
	  if(_scan[cryo_id][tpc_id].size()<=plane_id)
	    _scan[cryo_id][tpc_id].resize(plane_id+1);
	  _scan[cryo_id][tpc_id][plane_id] = _valid_nplanes;
	  ++_valid_nplanes;
	}
      }
    }
  }

  void SuperaWire::initialize()
  { SuperaBase::initialize(); }

  std::pair<size_t,size_t> SuperaWire::time_range(const geo::TPCGeo& tpc_geo,
						  const double x_min,
						  const double x_max)
  {
    LARCV_INFO() << "(xmin,xmax) = (" << x_min << "," << x_max << ")" << std::endl;
    auto pt = tpc_geo.ReferencePlane().GetCenter();
    pt.SetXYZ(x_min,pt.Y(),pt.Z());
    double dist_min = tpc_geo.DistanceFromReferencePlane(pt);
    pt.SetXYZ(x_max,pt.Y(),pt.Z());
    double dist_max = tpc_geo.DistanceFromReferencePlane(pt);
    if(dist_min > dist_max) std::swap(dist_min,dist_max);

    double min_time = (dist_min/supera::DriftVelocity() - supera::TriggerOffsetTPC()) / supera::TPCTickPeriod();
    double max_time = (dist_max/supera::DriftVelocity() - supera::TriggerOffsetTPC()) / supera::TPCTickPeriod();
    assert(min_time < max_time);

    LARCV_INFO() << "(tmin,tmax) = (" << min_time << "," << max_time << ")" << std::endl;

    double mid_time = min_time + (max_time - min_time) / 2.; 
    double num_time_half = _npx_rows * _time_compression / 2;

    size_t min_tick, max_tick;

    if( (mid_time + num_time_half) > supera::NumberTimeSamples() ) {
      max_tick = supera::NumberTimeSamples() - 1;
      min_tick = max_tick - _npx_rows * _time_compression;
      LARCV_INFO() << "mid_time + num_time_half (" << mid_time + num_time_half << ") abvoe max (" << supera::NumberTimeSamples() << ")"
		   << " ... tick range " << min_tick << " => " << max_tick <<std::endl;
    }
    else if( (mid_time + 1) > num_time_half ) {
      min_tick = 0;
      max_tick = _npx_rows * _time_compression - 1;
      LARCV_INFO() << "mid_time+1 (" << mid_time + 1 << ") too small"
		   << " ... tick range " << min_tick << " => " << max_tick <<std::endl;
    }else{
      max_tick = (int)(mid_time + num_time_half) - 1;
      min_tick = max_tick - _npx_rows + 1;
      LARCV_INFO() << " ... tick range " << min_tick << " => " << max_tick <<std::endl;
    }
    
    return std::pair<size_t,size_t>(min_tick,max_tick);

  }

  std::pair<size_t,size_t> SuperaWire::wire_range(const geo::PlaneGeo& plane_geo, 
						  const geo::Point_t& min_pt, 
						  const geo::Point_t& max_pt) 
  {
    LARCV_INFO() << "(ymin,ymax) = (" << min_pt.Y() << "," << max_pt.Y() << ")"
		 << " ... (zmin,zmax) = (" << min_pt.Z() << "," << max_pt.Z() << ")" << std::endl;

    static geo::Point_t pt0, pt1, pt2, pt3;
    auto const bbox = plane_geo.BoundingBox();
    /*
    pt0.SetXYZ(min_pt.X(),std::max(bbox.MinY(),min_pt.Y()),std::max(bbox.MinZ(),min_pt.Z()));
    pt1.SetXYZ(min_pt.X(),std::min(bbox.MaxY(),max_pt.Y()),std::max(bbox.MinZ(),min_pt.Z()));
    pt2.SetXYZ(max_pt.X(),std::max(bbox.MinY(),min_pt.Y()),std::min(bbox.MaxZ(),max_pt.Z()));
    pt3.SetXYZ(max_pt.X(),std::min(bbox.MaxY(),max_pt.Y()),std::min(bbox.MaxZ(),max_pt.Z()));
    auto wid0 = plane_geo.NearestWireID(pt0);
    auto wid1 = plane_geo.NearestWireID(pt1);
    auto wid2 = plane_geo.NearestWireID(pt2);
    auto wid3 = plane_geo.NearestWireID(pt3);
    size_t min_wire = std::min(wid0.Wire,wid1.Wire);
    min_wire = std::min(min_wire,(size_t)(wid2.Wire));
    min_wire = std::min(min_wire,(size_t)(wid3.Wire));

    size_t max_wire = std::max(wid0.Wire,wid1.Wire);
    max_wire = std::max(max_wire,(size_t)(wid2.Wire));
    max_wire = std::max(max_wire,(size_t)(wid3.Wire));
    */

    pt0.SetXYZ(bbox.CenterX(),min_pt.Y(),min_pt.Z());
    pt1.SetXYZ(bbox.CenterX(),max_pt.Y(),min_pt.Z());
    pt2.SetXYZ(bbox.CenterX(),min_pt.Y(),max_pt.Z());
    pt3.SetXYZ(bbox.CenterX(),max_pt.Y(),max_pt.Z());

    int wid0 = int(0.5 + plane_geo.WireCoordinate(pt0));
    int wid1 = int(0.5 + plane_geo.WireCoordinate(pt1));
    int wid2 = int(0.5 + plane_geo.WireCoordinate(pt2));
    int wid3 = int(0.5 + plane_geo.WireCoordinate(pt3));

    if(wid0 < 0) {wid0 = 0;} if(wid0 >= (int)(plane_geo.Nwires())) {wid0 = plane_geo.Nwires() - 1;}
    if(wid1 < 0) {wid1 = 0;} if(wid1 >= (int)(plane_geo.Nwires())) {wid1 = plane_geo.Nwires() - 1;}
    if(wid2 < 0) {wid2 = 0;} if(wid2 >= (int)(plane_geo.Nwires())) {wid2 = plane_geo.Nwires() - 1;}
    if(wid3 < 0) {wid3 = 0;} if(wid3 >= (int)(plane_geo.Nwires())) {wid3 = plane_geo.Nwires() - 1;}

    size_t min_wire = std::min(std::min(std::min(wid0,wid1),wid2),wid3);
    size_t max_wire = std::max(std::max(std::max(wid0,wid1),wid2),wid3);

    LARCV_INFO() << "(min_wire,max_wire) = (" << min_wire << "," << max_wire << ")" << std::endl;

    // Now define the range based on number of pixels requested + wire range
    size_t mid_wire = min_wire + (size_t)(((double)max_wire - (double)min_wire) / 2. + 0.5);

    size_t num_wire_half = _npx_columns / 2;
    if( (mid_wire + num_wire_half) > plane_geo.Nwires() ) {
      max_wire = plane_geo.Nwires()-1;
      min_wire = max_wire - _npx_columns + 1;
      LARCV_INFO() << "mid_wire + num_wire_half (" << mid_wire + num_wire_half << ") abvoe max (" << plane_geo.Nwires() << ")"
		   << " ... wire range " << min_wire << " => " << max_wire <<std::endl;
    }
    else if( (mid_wire+1) < num_wire_half ) {
      min_wire = 0;
      max_wire = _npx_columns -1;
      LARCV_INFO() << "mid_wire+1 (" << mid_wire + 1 << ") too small"
		   << " ... wire range " << min_wire << " => " << max_wire <<std::endl;
    }
    else{
      max_wire = (mid_wire + num_wire_half - 1);
      min_wire = mid_wire - num_wire_half;
      LARCV_INFO() << " ... wire range " << min_wire << " => " << max_wire <<std::endl;
    }
    return std::pair<size_t,size_t>(min_wire,max_wire);
  }


  bool SuperaWire::process(IOManager& mgr)
  {

    SuperaBase::process(mgr);
    larcv::Voxel3DMeta meta3d;

    if(!_ref_meta3d_cluster3d.empty()) {
      auto const& ev_cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(_ref_meta3d_cluster3d);
      meta3d = ev_cluster3d.meta();
    }
    else if(!_ref_meta3d_tensor3d.empty()) {
      auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_ref_meta3d_tensor3d);
      meta3d = ev_tensor3d.meta();
    }

    //
    // Define ImageMeta per plane, then store recob::Wire into VoxelSet per plane
    //
    // a) a list of voxel set to be filled
    std::vector<larcv::VoxelSet> vs_v(_valid_nplanes);
    // b) a list of meta2d 
    std::vector<larcv::ImageMeta> meta2d_v(_valid_nplanes);

    // First loop over configured planes and set the valid meta2d attributes (i.e. data)
    auto geop = lar::providerFrom<geo::Geometry>();
    for(unsigned int cryo_id=0; cryo_id<_scan.size(); ++cryo_id) {
      auto const& tpcs =  _scan.at(cryo_id);
      for(unsigned int tpc_id=0; tpc_id<tpcs.size(); ++tpc_id) {
	auto const& planes = tpcs.at(tpc_id);
	for(unsigned int plane_id=0; plane_id<planes.size(); ++plane_id) {
	  if(planes.at(plane_id)<0) continue;
	  LARCV_INFO() << "Creating meta: " << cryo_id << "-" << tpc_id << "_" << plane_id << std::endl;
	  // get the meta
	  auto& meta2d = meta2d_v.at(planes.at(plane_id));
	  // look up yz boundary
          auto const& plane_geo = geop->Plane({cryo_id, tpc_id, plane_id});
	  geo::Point_t min_pt3d, max_pt3d;
	  min_pt3d.SetXYZ(meta3d.bottom_left().x,meta3d.bottom_left().y,meta3d.bottom_left().z);
	  max_pt3d.SetXYZ(meta3d.top_right().x,meta3d.top_right().y,meta3d.top_right().z);
	  // Find the wire range
	  auto wrange = this->wire_range(plane_geo,min_pt3d,max_pt3d);
	  // Next look up x boundary
          auto const& tpc_geo = geop->TPC({cryo_id, tpc_id});
	  auto trange = this->time_range(tpc_geo,min_pt3d.X(),max_pt3d.X());
	  // Now you have x (wire) and y (time) range to store. Record in meta2d
	  meta2d = ImageMeta(wrange.first, trange.first, wrange.second+1, trange.second+1,
			     _npx_rows, _npx_columns, cryo_id*100+tpc_id*10+plane_id, larcv::kUnitCM);	  
	  LARCV_INFO() << meta2d.dump() << std::endl;
	}
      }
    }

    // Loop over wires to store actual image info
    auto const& wire_v = LArData<supera::LArWire_t>();
    
    for(auto const& wire : wire_v) {
      // Need ID info (cryo,tpc,wire num) to find relevant meta & voxelset
      auto ch  = wire.Channel();
      // To live with dumb detectors that use multiplex readout 
      auto wid_v = geop->ChannelToWire(ch);
      assert(wid_v.size() == 1);
      auto const& wid = wid_v[0];
      // Check to whether or not to store this
      auto const& idx = _scan.at(wid.Cryostat).at(wid.TPC).at(wid.Plane);
      if(idx < 0) continue;
      if(idx >= (int)(meta2d_v.size())) {
	LARCV_WARNING()<<idx << "Unexpected: " << meta2d_v.size() << " " << vs_v.size() << std::endl;
	throw larbys();
      }
      auto const& meta2d = meta2d_v.at(idx);
      if(wid.Wire < meta2d.min_x() || wid.Wire >= meta2d.max_x()) continue;
      auto& vs = vs_v.at(idx);
      
      // Loop over stored sparse vectors (can be multiple!)
      for (auto const& range : wire.SignalROI().get_ranges()) {
	
	auto const& adcs = range.data(); // actual pixel values for this wire (=column)
	int time_index = range.begin_index(); // start pixel time tick (=row)
	
	for(auto const& adc : adcs) {
	  if(adc <= _adc_threshold) continue;
	  if(time_index > meta2d.max_y()) break;
	  if(time_index < meta2d.min_y())
	    { time_index +=1; continue;}
	  vs.emplace(meta2d.id((double)wid.Wire,(double)(time_index)), adc, true);
	  time_index += 1;
	}
      }
    }

    // Store data
    for(unsigned int cryo_id=0; cryo_id<_scan.size(); ++cryo_id) {
      auto const& tpcs =  _scan.at(cryo_id);
      for(unsigned int tpc_id=0; tpc_id<tpcs.size(); ++tpc_id) {
	auto const& planes = tpcs.at(tpc_id);
	for(unsigned int plane_id=0; plane_id<planes.size(); ++plane_id) {
	  auto const& idx = planes.at(plane_id);
	  if(idx < 0) continue;
	  if(idx >= (int)(meta2d_v.size())) {
	    LARCV_WARNING()<<idx << "Unexpected: " << meta2d_v.size() << " " << vs_v.size() << std::endl;
	    throw larbys();
	  }
	  std::string output_name = _output_producer;
	  output_name = output_name + "_" + std::to_string(cryo_id);
	  output_name = output_name + "_" + std::to_string(tpc_id);
	  output_name = output_name + "_" + std::to_string(plane_id);
	  auto& output_event = mgr.get_data<larcv::EventSparseTensor2D>(output_name);
	  //output_event.emplace(std::move(vs_v.at(idx)),std::move(meta2d_v.at(idx)));
	  output_event.set(vs_v[idx],meta2d_v[idx]);
	  /*
	  std::cout<<cryo_id<<" "<<tpc_id<<" "<<plane_id 
		   <<" ... " << output_event.as_vector().front().size() 
		   << " " << output_event.as_vector().front().meta().id() << std::endl;
	  */
	}
      }
    }
    
    return true;
  }

  void SuperaWire::finalize()
  {}


}
#endif
