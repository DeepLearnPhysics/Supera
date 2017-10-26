/**
 * \file ParamsParticle.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class ParamsParticle
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAPARAMSParticle_H__
#define __SUPERAPARAMSParticle_H__

#include "FMWKInterface.h"

namespace supera {

  /**
     \class ParamsParticle
     User defined class ParamsParticle ... these comments are used to generate
     doxygen documentation!
  */
  class ParamsParticle {

  public:

    /// Default constructor
    ParamsParticle() {}

    /// Default destructor
    ~ParamsParticle() {}

    void configure(const supera::Config_t&);

    //
    // Getter
    //
    const std::string& OutParticleLabel()      const { return _out_roi_producer;      }

  private:

    std::string _out_roi_producer;

  };

}

#endif
/** @} */ // end of doxygen group

