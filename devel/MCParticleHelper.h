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

    ::larcv::Particle MakeParticle( const supera::LArMCParticle_t& mcp) const;

    ::larcv::Particle MakeParticle( const supera::LArMCParticle_t& mcp,
				    const larcv::Voxel3DMeta& meta3d) const;

  private:

    unsigned int _max_time_tick; ///< Maximum tick number in time
    unsigned int _time_padding;  ///< Padding in time axis (height) for MCParticleHelper::Format function
    unsigned int _wire_padding;  ///< Padding in wire axis (width) for MCParticleHelper::Format function
    bool _apply_sce;
  };
}

#endif
//#endif
//#endif
