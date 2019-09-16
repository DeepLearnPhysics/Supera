#ifndef __SUPERAMCPARTICLECLUSTER_DATA_H__
#define __SUPERAMCPARTICLECLUSTER_DATA_H__

#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace supera {

  enum ProcessType {
    kTrack,
    kNeutron,
    kPhoton,
    kPrimary,
    kCompton,
    kComptonHE,
    kDelta,
    kConversion,
    kIonization,
    kPhotoElectron,
    kDecay,
    kInvalidProcess
  };

  struct EDep {
    double x,y,z,t,e;
    EDep()
    { x = y = z = t = e = larcv::kINVALID_DOUBLE; }
  };
  
  struct EarlyPoints {
    std::vector<EDep> pts;
    void AddEDep(const EDep& pt) {
      if(pts.empty()) return;
      size_t loc=pts.size();
      while(1) {
	if(loc==0) break;
	size_t new_loc = loc-1;
	if(pts[new_loc].t < pt.t) break;
	loc = new_loc;
      }
      if(loc == pts.size()) return;
      for(int i=pts.size()-1; i>((int)(loc)); --i)
	pts[i] = pts[i-1];
      pts[loc] = pt;
      return;
    }
    const EDep& FirstPoint() const { return pts.front(); }
    void Direction(double& x, double& y, double& z) const 
    {
      double esum = 0;
      double dir[3] = {0., 0., 0.};
      if(pts.size()<2) {
	x = y = z = 0;
	return;
      }
      auto const& pt0 = pts.front();
      for(size_t i=1; i<pts.size(); ++i) {
	auto const& pt = pts[i];
	if(pt.x == larcv::kINVALID_DOUBLE) break;
	double norm = sqrt(pow(pt.x - pt0.x,2)+pow(pt.y - pt0.y,2)+pow(pt.z - pt0.z,2));
	dir[0] += ( pt.e * (pt.x - pt0.x) / norm);
	dir[1] += ( pt.e * (pt.y - pt0.y) / norm);
	dir[2] += ( pt.e * (pt.z - pt0.z) / norm);
	esum += pt.e;
      }
      x = dir[0] / esum;
      y = dir[1] / esum;
      z = dir[2] / esum;
    }
  };
  
  struct ParticleGroup {
    larcv::Particle part;
    std::vector<size_t> trackid_v;
    bool valid;
    ProcessType type;
    bool add_to_parent;
    larcv::VoxelSet vs;
    std::vector<larcv::VoxelSet> vs2d_v;
    EarlyPoints start;
    EDep last_pt;
    ParticleGroup(size_t num_planes=0)
    { valid=false; add_to_parent=false; type=kInvalidProcess; last_pt.t = -1.e9; start.pts.resize(10); vs2d_v.resize(num_planes);}
    void AddEDep(const EDep& pt)
    { if(pt.x == larcv::kINVALID_DOUBLE) return; start.AddEDep(pt); if(pt.t > last_pt.t) last_pt = pt; }

    void Merge(ParticleGroup& child) {
      for(auto const& vox : child.vs.as_vector())
	this->vs.emplace(vox.id(),vox.value(),true);
      for(auto const& pt : child.start.pts)
	this->AddEDep(pt);
      this->AddEDep(child.last_pt);
      this->trackid_v.push_back(child.part.track_id());
      for(auto const& trackid : child.trackid_v)
	this->trackid_v.push_back(trackid);
      for(size_t plane_id=0; plane_id < vs2d_v.size(); ++plane_id) {
	auto& vs2d = vs2d_v[plane_id];
	auto& child_vs2d = child.vs2d_v[plane_id];
	for(auto const& vox : child_vs2d.as_vector())
	  vs2d.emplace(vox.id(),vox.value(),true);
	child_vs2d.clear_data();
      }
      child.vs.clear_data();
      child.valid=false;
    }
  };
}

#endif
