/**
 * \file SuperaWire.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaWire
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAWIRE_H__
#define __SUPERAWIRE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "SuperaBase.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaWire ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaWire : public SuperaBase {

  public:
    /// Default constructor
    SuperaWire(const std::string name = "SuperaWire");

    /// Default destructor
    ~SuperaWire() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    int plane_index(unsigned int cryo_id, unsigned int tpc_id, unsigned int plane_id);
    std::pair<size_t, size_t> time_range(const geo::TPCGeo& tpc_geo,
                                         const double x_min,
                                         const double x_max);
    std::pair<size_t, size_t> wire_range(const geo::PlaneGeo& plane_geo,
                                         const geo::Point_t& min_pt,
                                         const geo::Point_t& max_pt);
    std::vector<std::vector<std::vector<int>>> _scan;
    std::string _output_producer;
    size_t _valid_nplanes;
    std::string _ref_meta3d_cluster3d;
    std::string _ref_meta3d_tensor3d;
    int _npx_rows, _npx_columns;
    double _time_compression;
    double _adc_threshold;
  };

  /**
     \class larcv::SuperaWireFactory
     \brief A concrete factory class for larcv::SuperaWire
  */
  class SuperaWireProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaWireProcessFactory() { ProcessFactory::get().add_factory("SuperaWire", this); }
    /// dtor
    ~SuperaWireProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaWire(instance_name); }
  };

}
#endif
//#endif
//#endif
/** @} */ // end of doxygen group
