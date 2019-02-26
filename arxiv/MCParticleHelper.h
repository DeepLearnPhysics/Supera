#ifndef __SUPERA_MCPARTICLEHELPER_H__
#define __SUPERA_MCPARTICLEHELPER_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include <vector>

// LArSoft
//#include "MCBase/MCTrack.h"
//#include "MCBase/MCShower.h"
//#include "Simulation/SimChannel.h"
//#include "SimulationBase/MCParticle.h"

// LArCV
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Particle.h"
#include "larcv/core/DataFormat/BBox.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
#include "FMWKInterface.h"
#include "Range.h"
namespace supera {

  typedef ::larcv::Range<unsigned int> WTRange_t;
  typedef std::vector<supera::WTRange_t> WTRangeArray_t;

  class MCParticleHelper : public ::larcv::larcv_base {

  public:

    MCParticleHelper() : larcv::larcv_base("MCParticleHelper")
      , _max_time_tick(9600)
      , _time_padding(10)
      , _wire_padding(10)
    {}

    virtual ~MCParticleHelper() {}

    void configure(const supera::Config_t& cfg);
    /**
       Given single MCTrack, returns length 4 range array (3 planes + time) \n
       which contains all trajectory points of input MCTrack.
    */
    WTRangeArray_t WireTimeBoundary( const supera::LArMCTrack_t& mct ) const;

    WTRangeArray_t WireTimeBoundary( const supera::LArMCTrack_t& mct,
                                     const std::vector<supera::LArSimCh_t>& sch_v ) const;

    WTRangeArray_t WireTimeBoundary( const supera::LArMCShower_t& mcs ) const;

    WTRangeArray_t WireTimeBoundary( const supera::LArMCShower_t& mcs,
                                     const std::vector<supera::LArSimCh_t>& sch_v ) const;

    std::vector<larcv::BBox2D> MakeBBox2D( const supera::LArMCTrack_t& mct,
                                           const int time_offset ) const;

    std::vector<larcv::BBox2D> MakeBBox2D( const supera::LArMCTrack_t& mct,
                                           const std::vector<supera::LArSimCh_t>& sch_v,
                                           const int time_offset ) const;

    std::vector<larcv::BBox2D> MakeBBox2D( const supera::LArMCShower_t& mcs,
                                           const int time_offset ) const;

    std::vector<larcv::BBox2D> MakeBBox2D( const supera::LArMCShower_t& mcs,
                                           const std::vector<supera::LArSimCh_t>& sch_v,
                                           const int time_offset ) const;

    ::larcv::Particle MakeParticle( const supera::LArMCTrack_t& mct) const;

    ::larcv::Particle MakeParticle( const supera::LArMCTrack_t& mct,
				    const larcv::Voxel3DMeta& meta3d) const;

    ::larcv::Particle MakeParticle( const supera::LArMCShower_t& mcs) const;

    ::larcv::Particle MakeParticle( const supera::LArMCTrack_t& mct,
				    const std::vector<supera::LArSimCh_t>& sch_v,
				    const int time_offset) const;

    ::larcv::Particle MakeParticle( const supera::LArMCShower_t& mcs,
				    const std::vector<supera::LArSimCh_t>& sch_v,
				    const int time_offset ) const;

    ::larcv::Particle MakeParticle( const supera::LArMCTrack_t& mct,
				    const std::vector<supera::LArSimCh_t>& sch_v,
				    const int time_offset,
				    const larcv::Voxel3DMeta& meta3d) const;

    ::larcv::Particle MakeParticle( const supera::LArMCShower_t& mcs,
				    const std::vector<supera::LArSimCh_t>& sch_v,
				    const int time_offset,
				    const larcv::Voxel3DMeta& meta3d) const;

  private:

    std::vector<larcv::BBox2D> WTRange2BB(const WTRangeArray_t&) const;

    unsigned int _max_time_tick; ///< Maximum tick number in time
    unsigned int _time_padding;  ///< Padding in time axis (height) for MCParticleHelper::Format function
    unsigned int _wire_padding;  ///< Padding in wire axis (width) for MCParticleHelper::Format function
    bool _apply_sce;
  };
}

#endif
//#endif
//#endif
