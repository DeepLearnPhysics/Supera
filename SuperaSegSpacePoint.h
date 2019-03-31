/**
 * \file SuperaSegSpacePoint.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaSegSpacePoint
 *
 * @author laura
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERASEGSPACEPOINT_H__
#define __SUPERASEGSPACEPOINT_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaSegSpacePoint ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaSegSpacePoint : public SuperaBase {

  public:

    /// Default constructor
    SuperaSegSpacePoint(const std::string name = "SuperaSegSpacePoint");

    /// Default destructor
    ~SuperaSegSpacePoint() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    std::string _output_label, _data_label;
		float _distance_threshold;
  };

  /**
     \class larcv::SuperaSegSpacePointFactory
     \brief A concrete factory class for larcv::SuperaSegSpacePoint
  */
  class SuperaSegSpacePointProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaSegSpacePointProcessFactory() { ProcessFactory::get().add_factory("SuperaSegSpacePoint", this); }
    /// dtor
    ~SuperaSegSpacePointProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaSegSpacePoint(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group

