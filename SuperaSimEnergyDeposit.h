/**
 * \file SuperaSimEnergyDeposit.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaSimEnergyDeposit
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERASIMENERGYDEPOSIT_H__
#define __SUPERASIMENERGYDEPOSIT_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"
#include "larcv/core/DataFormat/BBox.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaSimEnergyDeposit ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaSimEnergyDeposit : public SuperaBase {

  public:

    /// Default constructor
    SuperaSimEnergyDeposit(const std::string name = "SuperaSimEnergyDeposit");

    /// Default destructor
    ~SuperaSimEnergyDeposit() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    bool _store_dx, _store_dq, _store_dp, _store_dt, _store_at, _store_dedx;
    std::string _output_label, _particle_label;
    BBox3D _world_bounds;

  };

  /**
     \class larcv::SuperaSimEnergyDepositFactory
     \brief A concrete factory class for larcv::SuperaSimEnergyDeposit
  */
  class SuperaSimEnergyDepositProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaSimEnergyDepositProcessFactory() { ProcessFactory::get().add_factory("SuperaSimEnergyDeposit", this); }
    /// dtor
    ~SuperaSimEnergyDepositProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaSimEnergyDeposit(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group

