/**
 * \file SuperaToy.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaToy
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERATOY_H__
#define __SUPERATOY_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaToy ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaToy : public SuperaBase {

  public:

    /// Default constructor
    SuperaToy(const std::string name = "SuperaToy");

    /// Default destructor
    ~SuperaToy() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();
    
  private:

    larcv::Voxel3DMeta get_meta3d(IOManager& mgr) const;
    double _nsigma_match_time;
    bool _debug, _use_true_pos;
    std::string _hit_producer, _sps_producer, _ref_meta3d_cluster3d, _ref_meta3d_tensor3d;
  };

  /**
     \class larcv::SuperaToyFactory
     \brief A concrete factory class for larcv::SuperaToy
  */
  class SuperaToyProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaToyProcessFactory() { ProcessFactory::get().add_factory("SuperaToy", this); }
    /// dtor
    ~SuperaToyProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaToy(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group

