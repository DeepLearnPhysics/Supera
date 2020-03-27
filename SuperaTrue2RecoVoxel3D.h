/**
 * \file SuperaTrue2RecoVoxel3D.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaTrue2RecoVoxel3D
 *
 * @author Laura Domine, Patrick Tsang
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SuperaTrue2RecoVoxel3D_H__
#define __SuperaTrue2RecoVoxel3D_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include <functional>
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaTrue2RecoVoxel3D ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaTrue2RecoVoxel3D : public SuperaBase {

  public:

    /// Default constructor
    SuperaTrue2RecoVoxel3D(const std::string name = "SuperaTrue2RecoVoxel3D");

    /// Default destructor
    ~SuperaTrue2RecoVoxel3D() {}

    void configure(const PSet&);

    void initialize();

	 	bool process(IOManager& mgr);

    void finalize();

  private:

    larcv::Voxel3DMeta get_meta3d(IOManager& mgr) const;
    std::string _hit_producer, _sps_producer, _ref_meta3d_cluster3d, _ref_meta3d_tensor3d;
    std::string _output_tensor3d, _output_cluster3d;
    bool _debug, _use_true_pos;
    bool _twofold_matching;
    bool _dump_to_csv;

    bool _post_averaging;
    double _post_averaging_threshold;

    bool _hit_peak_finding;
    double _hit_threshold_ne;
    double _hit_window_ticks;

    // composite index of (voxel_id, track_id)
    struct TrackVoxel_t
    {
      unsigned int voxel_id;
      int track_id;

      TrackVoxel_t(unsigned int v, int t) : voxel_id(v), track_id(t) {}
      friend inline bool operator< (const TrackVoxel_t &lhs, const TrackVoxel_t &rhs){
        if (lhs.voxel_id == rhs.voxel_id)
          return lhs.track_id < rhs.track_id;
        return lhs.voxel_id < rhs.voxel_id;
      }
    };

    // minimal hit info (time, voxel_id, track_id)
    struct TrueHit_t
    {
      double time;
      friend inline bool operator<(const TrueHit_t &lhs, const TrueHit_t &rhs){
        return lhs.time < rhs.time;
      }
      std::vector<TrackVoxel_t> track_voxel_ids;
      std::vector<double> n_electrons;
    };

    // find true hits in [t_start, t_end]
    void find_hits(const std::vector<TrueHit_t>& hits, double t_start, double t_end,
        std::set<TrackVoxel_t>& track_voxel_ids);

    // find peak (per track_id) of true hits in [t_start, t_end]
    void find_hit_peaks(const std::vector<TrueHit_t>& hits, double t_start, double t_end,
        std::set<TrackVoxel_t>& track_voxel_ids);

    // make ghost labels with simple overlapping 2 or 3 true hits
    void set_ghost(
        const std::map<unsigned int, std::set<TrackVoxel_t>>& reco2true,
        std::map<unsigned int, bool>& ghosts,
        std::function<void(unsigned int, unsigned int)> const& insert_true2reco);

    // make ghost labels with averaging over reco pts for each true pt
    void set_ghost_with_averaging(
        const std::map<unsigned int, std::set<TrackVoxel_t>>& reco2true,
        const larcv::Voxel3DMeta& meta3d,
        std::map<unsigned int, bool>& ghosts,
        std::function<void(unsigned int, unsigned int)> const& insert_true2reco);
  };

  /**
     \class larcv::SuperaTrue2RecoVoxel3DFactory
     \brief A concrete factory class for larcv::SuperaTrue2RecoVoxel3D
  */
  class SuperaTrue2RecoVoxel3DProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaTrue2RecoVoxel3DProcessFactory() { ProcessFactory::get().add_factory("SuperaTrue2RecoVoxel3D", this); }
    /// dtor
    ~SuperaTrue2RecoVoxel3DProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaTrue2RecoVoxel3D(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
