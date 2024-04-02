/**
 * \file SuperaOptical.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaOptical
 *
 * @author Temigo
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAOPTICAL_H__
#define __SUPERAOPTICAL_H__
#include "SuperaBase.h"
#include <vector>

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaOptical ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaOptical : public SuperaBase {

  public:
    /// Default constructor
    SuperaOptical(const std::string name = "SuperaOptical");

    /// Default destructor
    ~SuperaOptical() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::vector<std::string> _opflash_producer_label_v;
    std::vector<std::string> _opflash_output_label_v;
  };

  /**
     \class larcv::SuperaOpticalFactory
     \brief A concrete factory class for larcv::SuperaOptical
  */
  class SuperaOpticalProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaOpticalProcessFactory() { ProcessFactory::get().add_factory("SuperaOptical", this); }
    /// dtor
    ~SuperaOpticalProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaOptical(instance_name);
    }
  };

}
#endif
/** @} */ // end of doxygen group
