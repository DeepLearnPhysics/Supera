/**
 * \file SuperaMCTruth.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaMCTruth
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMCTRUTH_H__
#define __SUPERAMCTRUTH_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "SuperaBase.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaMCTruth ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaMCTruth : public SuperaBase {

  public:
    /// Default constructor
    SuperaMCTruth(const std::string name = "SuperaMCTruth");

    /// Default destructor
    ~SuperaMCTruth() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    unsigned short _pass_origin;
    std::string _output_label;
    std::vector<std::string> _producer_labels;
  };

  /**
     \class larcv::SuperaMCTruthFactory
     \brief A concrete factory class for larcv::SuperaMCTruth
  */
  class SuperaMCTruthProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCTruthProcessFactory() { ProcessFactory::get().add_factory("SuperaMCTruth", this); }
    /// dtor
    ~SuperaMCTruthProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaMCTruth(instance_name);
    }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
