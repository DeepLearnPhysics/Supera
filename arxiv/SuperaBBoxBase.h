/**
 * \file SuperaBBoxBase.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaBBoxBase
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERABBOXBASE_H__
#define __SUPERABBOXBASE_H__
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "ImageMetaMakerBase.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaBBoxBase ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaBBoxBase : public SuperaBase {

  public:
    
    /// Default constructor
    SuperaBBoxBase(const std::string name="SuperaBBoxBase");
    
    /// Default destructor
    ~SuperaBBoxBase(){ if(_meta_maker) delete _meta_maker; }

    bool is(const std::string question) const;

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    supera::ImageMetaMakerBase* _meta_maker;
    std::map<supera::RSEID,std::array<double,3> > _constraint_m;
  };

  /**
     \class larcv::SuperaBBoxBaseFactory
     \brief A concrete factory class for larcv::SuperaBBoxBase
  */
  class SuperaBBoxBaseProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaBBoxBaseProcessFactory() { ProcessFactory::get().add_factory("SuperaBBoxBase",this); }
    /// dtor
    ~SuperaBBoxBaseProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaBBoxBase(instance_name); }
  };

}
#endif
/** @} */ // end of doxygen group 

