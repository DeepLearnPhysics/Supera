/**
 * \file SuperaSpacePoint.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaSpacePoint
 *
 * @author kvtsang
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERASPACEPOINT_H__
#define __SUPERASPACEPOINT_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"
#include "larcv/core/DataFormat/BBox.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaSpacePoint ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaSpacePoint : public SuperaBase {

  public:
    /// Default constructor
    SuperaSpacePoint(const std::string name = "SuperaSpacePoint");

    /// Default destructor
    ~SuperaSpacePoint() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    static float get_common_charge(const std::vector<art::Ptr<recob::Hit>>& hits);

  private:
    larcv::BBox3D _world_bounds;
    std::vector<std::string> _producer_labels;
    std::string _output_label;
    size_t _max_debug_dropping = 0; // Max debug message for dropping space points
    unsigned short _n_planes = 3;
    float _shift_x = 0.;
    bool _store_wire_info = false;
    std::unordered_set<std::string> _drop_output;
    std::vector<double> _reco_charge_range;
  };

  /**
     \class larcv::SuperaSpacePointFactory
     \brief A concrete factory class for larcv::SuperaSpacePoint
  */
  class SuperaSpacePointProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaSpacePointProcessFactory()
    {
      ProcessFactory::get().add_factory("SuperaSpacePoint", this);
    }
    /// dtor
    ~SuperaSpacePointProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaSpacePoint(instance_name);
    }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
