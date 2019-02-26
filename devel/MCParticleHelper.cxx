#ifndef __SUPERA_MCPARTICLEHELPER_CXX__
#define __SUPERA_MCPARTICLEHELPER_CXX__

#include "MCParticleHelper.h"
#include "larcv/core/Base/larbys.h"
#include "FMWKInterface.h"
#include <TLorentzVector.h> // ROOT
#include <set>
namespace supera {

  void MCParticleHelper::configure(const supera::Config_t& cfg)
  {
    LARCV_DEBUG() << "start" << std::endl;
    set_verbosity((::larcv::msg::Level_t)(cfg.get<unsigned short>("Verbosity", logger().level())));
    _max_time_tick = cfg.get<unsigned int>("MaxTimeTick");
    _time_padding  = cfg.get<unsigned int>("TimePadding");
    _wire_padding  = cfg.get<unsigned int>("WirePadding");
    _apply_sce     = cfg.get<bool>("ApplySCE");

    LARCV_NORMAL() << "Configuration called..." << std::endl
                   << " Padding   (wire,time) = (" << _wire_padding << "," << _time_padding << ")" << std::endl;
  }

  std::vector<larcv::Particle> MCParticleHelper::Convert( const std::vector<supera::LArMCParticle_t>& mcp_v) const
  {
    LARCV_DEBUG() << "start" << std::endl;

    std::vector<larcv::Particle> res_v;
    res_v.reserve(mcp_v.size());
    for(size_t idx=0; idx<mcp_v.size(); ++idx) {
      auto const& mcp = mcp_v.at(idx);
      
      LARCV_INFO() << "Assessing MCTrack G4Track ID = " << mcp.TrackId() << " PdgCode " << mcp.PdgCode() << std::endl;

      larcv::Particle res;
      double xyz[3] = {0.};
      
      
      
      res.mcst_index ( idx           ); 
      res.track_id   ( mcp.TrackId() );
      res.pdg_code   ( mcp.PdgCode() );
      res.momentum   ( 


    inline void id              (InstanceID_t id  )  { _id = id;         }
    inline void mcst_index      (MCSTIndex_t id )    { _mcst_index = id;    }
    inline void mct_index       (MCTIndex_t id )     { _mct_index = id;     }
    inline void shape           (ShapeType_t shape ) { _shape = shape;      }
    inline void nu_current_type (short curr) {_current_type = curr; }
    inline void nu_interaction_type (short itype) {_interaction_type = itype; }
    // particle's info setter                                                                                                              
    inline void track_id        (unsigned int id )   { _trackid = id;       }
    inline void pdg_code        (int code)           { _pdg = code;         }
    inline void momentum        (double px, double py, double pz) { _px = px; _py = py; _pz = pz; }
    inline void position        (const larcv::Vertex& vtx) { _vtx = vtx;    }
    inline void position        (double x, double y, double z, double t) { _vtx = Vertex(x,y,z,t); }
    inline void end_position    (const larcv::Vertex& vtx) { _end_pt = vtx; }
    inline void end_position    (double x, double y, double z, double t) { _end_pt = Vertex(x,y,z,t); }
    inline void first_step      (const larcv::Vertex& vtx) { _first_step = vtx; }
    inline void first_step      (double x, double y, double z, double t) { _first_step = Vertex(x,y,z,t); }
    inline void last_step       (const larcv::Vertex& vtx) { _last_step = vtx; }
    inline void last_step       (double x, double y, double z, double t) { _last_step = Vertex(x,y,z,t); }
    inline void distance_travel ( double dist ) { _dist_travel = dist; }
    inline void energy_init     (double e)           { _energy_init = e;    }
    inline void energy_deposit  (double e)           { _energy_deposit = e; }
    inline void creation_process (const std::string& proc) { _process = proc; }
    inline void boundingbox_2d(const std::vector<larcv::BBox2D>& bb_v) { _bb2d_v = bb_v; }
    inline void boundingbox_2d(const BBox2D& bb, ProjectionID_t id) { _bb2d_v.resize(id+1); _bb2d_v[id] = bb; }
    inline void boundingbox_3d(const BBox3D& bb) { _bb3d = bb; }
    inline void num_voxels(int count) { _num_voxels = count; }
    //inline void type_score (const std::vector<float>& score_v) { _type_score_v = score_v; }                                              
    // parent info setter                                                                                                                  
    inline void parent_track_id (unsigned int id )   { _parent_trackid = id;}
    inline void parent_pdg_code (int code)           { _parent_pdg = code;  }
    inline void parent_position (const larcv::Vertex& vtx) { _parent_vtx = vtx; }
    inline void parent_position (double x, double y, double z, double t) { _parent_vtx = Vertex(x,y,z,t); }
    // ancestor info setter                                                                                                                
    inline void ancestor_track_id (unsigned int id )   { _ancestor_trackid = id;}
    inline void ancestor_pdg_code (int code)           { _ancestor_pdg = code;  }
    inline void ancestor_position (const larcv::Vertex& vtx) { _ancestor_vtx = vtx; }
    inline void ancestor_position (double x, double y, double z, double t) { _ancestor_vtx = Vertex(x,y,z,t); }











    ::larcv::Particle res;
    res.shape(::larcv::kShapeTrack);
    //res.Type(::larcv::PdgCode2ROIType(mct.PdgCode()));
    if (mct.size())
      res.energy_deposit(mct.front().E() - mct.back().E());
    else
      res.energy_deposit(0);
    res.energy_init(mct.Start().E());

    xyz[0] = mct.Start().X();
    xyz[1] = mct.Start().Y();
    xyz[2] = mct.Start().Z();
    if (_apply_sce) ApplySCE(xyz);
    res.position(xyz[0], xyz[1], xyz[2], mct.Start().T());

    xyz[0] = mct.End().X();
    xyz[1] = mct.End().Y();
    xyz[2] = mct.End().Z();
    if (_apply_sce) ApplySCE(xyz);
    res.end_position(xyz[0], xyz[1], xyz[2], mct.End().T());

    res.creation_process(mct.Process());

    if (mct.size() > 0) {
      auto const& first_step = mct.front();
      xyz[0] = first_step.X();
      xyz[1] = first_step.Y();
      xyz[2] = first_step.Z();
      if (_apply_sce) ApplySCE(xyz);
      res.first_step(xyz[0], xyz[1], xyz[2], first_step.T());
    }
    if (mct.size() > 1) {
      auto const& last_step = mct.back();
      res.last_step(last_step.X(), last_step.Y(), last_step.Z(), last_step.T());
      double length = 0;
      for (size_t step_idx = 1; step_idx < mct.size(); ++step_idx) {
        auto const& step1 = mct[step_idx - 1];
        auto const& step2 = mct[step_idx];
        length += sqrt(pow(step1.X() - step2.X(), 2) + pow(step1.Y() - step2.Y(), 2) + pow(step1.Z() - step2.Z(), 2));
      }
      res.distance_travel(length);
    }
    res.momentum(mct.Start().Px(), mct.Start().Py(), mct.Start().Pz());
    res.pdg_code(mct.PdgCode());
    res.parent_pdg_code(mct.MotherPdgCode());
    res.track_id(mct.TrackID());
    res.parent_track_id(mct.MotherTrackID());

    xyz[0] = mct.MotherStart().X();
    xyz[1] = mct.MotherStart().Y();
    xyz[2] = mct.MotherStart().Z();
    if (_apply_sce) ApplySCE(xyz);

    res.parent_position(xyz[0], xyz[1], xyz[2], mct.MotherStart().T());
    /*
    res.parent_momentum(mct.MotherStart().Px(),
                        mct.MotherStart().Py(),
                        mct.MotherStart().Pz());
    */
    LARCV_INFO() << res.dump();
    return res;
  }

}
#endif

// Local Variables:
// mode: c++
// End:
