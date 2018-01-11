/**
 * \file SuperaInstanceImage.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaInstanceImage
 *
 * @author taritree
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERA_INSTANCE_IMAGE_H__
#define __SUPERA_INSTANCE_IMAGE_H__
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "ImageMetaMaker.h"
#include "ParamsImage2D.h"
#include "larcv/core/DataFormat/Image2D.h"
#include <string>

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaInstanceImage ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaInstanceImage : public SuperaBase,
    public supera::ParamsImage2D,
    public supera::ImageMetaMaker {
    
  public:
    
    /// Default constructor
    SuperaInstanceImage(const std::string name="SuperaInstanceImage");
    
    /// Default destructor
    ~SuperaInstanceImage(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    unsigned short _origin;
    std::string m_ancestor_label;
  };

  /**
     \class larcv::SuperaInstanceImageFactory
     \brief A concrete factory class for larcv::SuperaInstanceImage
  */
  class SuperaInstanceImageProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaInstanceImageProcessFactory() { ProcessFactory::get().add_factory("SuperaInstanceImage",this); }
    /// dtor
    ~SuperaInstanceImageProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaInstanceImage(instance_name); }
  };

}
#endif
/** @} */ // end of doxygen group 

