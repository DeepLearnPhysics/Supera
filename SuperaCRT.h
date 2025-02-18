/**
 * \file SuperaCRT.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaCRT
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
     User defined class SuperaCRT ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaCRT : public SuperaBase {

  public:
    /// Default constructor
    SuperaCRT(const std::string name = "SuperaCRT");

    /// Default destructor
    ~SuperaCRT() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::vector<std::string> _crthit_producer_label_v;
    std::vector<std::string> _crthit_output_label_v;
  };

  /**
     \class larcv::SuperaCRTFactory
     \brief A concrete factory class for larcv::SuperaCRT
  */
  class SuperaCRTProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaCRTProcessFactory() { ProcessFactory::get().add_factory("SuperaCRT", this); }
    /// dtor
    ~SuperaCRTProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaCRT(instance_name); }
  };

}

#endif
