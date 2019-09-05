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
    EarlyPoints start;
    EDep last_pt;
    ParticleGroup()
    { valid=false; add_to_parent=false; type=kInvalidProcess; last_pt.t = -1.e9; start.pts.resize(10);}
    void AddEDep(const EDep& pt)
    { if(pt.x == larcv::kINVALID_DOUBLE) return; start.AddEDep(pt); if(pt.t > last_pt.t) last_pt = pt; }
  };
}

#endif
