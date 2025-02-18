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
#include "FMWKInterface.h"
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Particle.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
namespace supera {

  class MCParticleHelper : public ::larcv::larcv_base {

  public:
    MCParticleHelper() : larcv::larcv_base("MCParticleHelper") {}

    virtual ~MCParticleHelper() {}

    void configure(const supera::Config_t& cfg);
    ::larcv::Particle MakeParticle(const supera::LArMCTrack_t& mct,
                                   const larcv::Voxel3DMeta& meta) const;

    ::larcv::Particle MakeParticle(const supera::LArMCShower_t& mcs) const;

  private:
    bool _apply_sce;
  };
}

#endif
//#endif
//#endif
