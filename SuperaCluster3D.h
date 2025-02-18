/**
 * \file SuperaCluster3D.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaCluster3D
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERACLUSTER3D_H__
#define __SUPERACLUSTER3D_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaCluster3D ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaCluster3D : public SuperaBase {

  public:
    /// Default constructor
    SuperaCluster3D(const std::string name = "SuperaCluster3D");

    /// Default destructor
    ~SuperaCluster3D() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::string _output_label;
  };

  /**
     \class larcv::SuperaCluster3DFactory
     \brief A concrete factory class for larcv::SuperaCluster3D
  */
  class SuperaCluster3DProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaCluster3DProcessFactory() { ProcessFactory::get().add_factory("SuperaCluster3D", this); }
    /// dtor
    ~SuperaCluster3DProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaCluster3D(instance_name);
    }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
