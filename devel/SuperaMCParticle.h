/**
 * \file SuperaMCParticle.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaMCParticle
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMCPARTICLE_H__
#define __SUPERAMCPARTICLE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaMCParticle ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaMCParticle : public SuperaBase {

  public:

    /// Default constructor
    SuperaMCParticle(const std::string name = "SuperaMCParticle");

    /// Default destructor
    ~SuperaMCParticle() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    std::string _cluster3d_labels;
    std::string _tensor3d_labels;
  };

  /**
     \class larcv::SuperaMCParticleFactory
     \brief A concrete factory class for larcv::SuperaMCParticle
  */
  class SuperaMCParticleProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCParticleProcessFactory() { ProcessFactory::get().add_factory("SuperaMCParticle", this); }
    /// dtor
    ~SuperaMCParticleProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaMCParticle(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group

