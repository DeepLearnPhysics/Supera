/**
 * \file SegSpacePoint.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SegSpacePoint
 *
 * @author laura
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SEGSPACEPOINT_H__
#define __SEGSPACEPOINT_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SegSpacePoint ... these comments are used to generate
     doxygen documentation!
  */
  class SegSpacePoint : public ProcessBase {

  public:
    /// Default constructor
    SegSpacePoint(const std::string name = "SegSpacePoint");

    /// Default destructor
    ~SegSpacePoint() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::string _output_label, _data_label;
    float _distance_threshold;
  };

  /**
     \class larcv::SegSpacePointFactory
     \brief A concrete factory class for larcv::SegSpacePoint
  */
  class SegSpacePointProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SegSpacePointProcessFactory() { ProcessFactory::get().add_factory("SegSpacePoint", this); }
    /// dtor
    ~SegSpacePointProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SegSpacePoint(instance_name);
    }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
