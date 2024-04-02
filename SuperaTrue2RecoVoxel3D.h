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
#include "FMWKInterface.h"
#include "RecoVoxel3D.h"
#include "SuperaBase.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
#include <functional>

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

    const std::unordered_set<TrackVoxel_t>& find_true(VoxelID_t reco_id) const;

    const std::unordered_set<RecoVoxel3D>& find_reco(int track_id, VoxelID_t true_id) const;

    const auto& get_reco2true() const { return _reco2true; };
    const auto& get_true2reco() const { return _true2reco; };
    const auto get_true2reco_bytrack() const { return this->contract_true2reco_bytrack(); }

  private:
    larcv::Voxel3DMeta get_meta3d(IOManager& mgr) const;
    std::string _ref_meta3d_cluster3d, _ref_meta3d_tensor3d;
    std::vector<std::string> _sps_producer_v;
    std::string _output_tensor3d, _output_cluster3d;
    bool _debug, _use_true_pos;
    bool _twofold_matching;
    bool _dump_to_csv;
    bool _useOrigTrackID; ///< Whether to use origTrackID or trackID from SimChannel

    bool _post_averaging;
    double _post_averaging_threshold;

    bool _hit_peak_finding;
    double _hit_threshold_ne;
    double _hit_window_ticks;
    double _voxel_size_factor;
    double _voxel_distance_threshold;

    std::vector<double> _reco_charge_range;

    // Map (true_id, track_id) -> [RecoVoxel3D]
    std::map<TrackVoxel_t, std::unordered_set<RecoVoxel3D>> _true2reco;

    // reverse map for RecoVoxel3D -> [true_id, track_id]
    std::map<RecoVoxel3D, std::unordered_set<TrackVoxel_t>> _reco2true;

    // empty set for lookup return with invalid key
    std::unordered_set<TrackVoxel_t> _empty_true;
    std::unordered_set<RecoVoxel3D> _empty_reco;

    // helper function to build reco2true and true2reco maps
    template <typename K, typename V>
    inline void insert_one_to_many(std::map<K, std::unordered_set<V>>& m,
                                   K const& key,
                                   V const& value) const
    {
      auto itr = m.find(key);
      if (itr == m.end())
        m.emplace(key, std::unordered_set<V>({value}));
      else
        itr->second.insert(value);
    }

    // clear contents of reco2true and true2reco maps
    void clear_maps();

    // find true hits in [t_start, t_end]
    void find_hits(const std::vector<TrueHit_t>& hits,
                   double t_start,
                   double t_end,
                   std::set<TrackVoxel_t>& track_voxel_ids);

    // find peak (per track_id) of true hits in [t_start, t_end]
    void find_hit_peaks(const std::vector<TrueHit_t>& hits,
                        double t_start,
                        double t_end,
                        std::set<TrackVoxel_t>& track_voxel_ids);

    // make ghost labels with simple overlapping 2 or 3 true hits
    void set_ghost();

    // make ghost labels with averaging over reco pts for each true pt
    void set_ghost_with_averaging(const larcv::Voxel3DMeta& meta3d);

    void set_intersection_distance(const larcv::Voxel3DMeta& meta3d,
                                   const std::set<TrackVoxel_t>& v1,
                                   const std::set<TrackVoxel_t>& v2,
                                   std::set<TrackVoxel_t>& overlaps);

    void set_intersection_factor(const larcv::Voxel3DMeta& meta3d,
                                 const std::set<TrackVoxel_t>& v1,
                                 const std::set<TrackVoxel_t>& v2,
                                 std::set<TrackVoxel_t>& overlaps);

    bool is_ghost(VoxelID_t reco_id) const;

    // build a map of reco_id <-> true_id, ignoring track_id
    std::map<VoxelID_t, std::unordered_set<VoxelID_t>> contract_true2reco();
    std::map<VoxelID_t, std::unordered_set<VoxelID_t>> contract_reco2true();
    std::vector<std::map<VoxelID_t, std::unordered_set<VoxelID_t>>> contract_true2reco_bytrack()
      const;
  };

  /**
     \class larcv::SuperaTrue2RecoVoxel3DFactory
     \brief A concrete factory class for larcv::SuperaTrue2RecoVoxel3D
  */
  class SuperaTrue2RecoVoxel3DProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaTrue2RecoVoxel3DProcessFactory()
    {
      ProcessFactory::get().add_factory("SuperaTrue2RecoVoxel3D", this);
    }
    /// dtor
    ~SuperaTrue2RecoVoxel3DProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaTrue2RecoVoxel3D(instance_name);
    }
  };
}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
