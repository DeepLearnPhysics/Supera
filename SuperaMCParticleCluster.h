/**
 * \file SuperaMCParticleCluster.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaMCParticleCluster
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMCPARTICLECLUSTER_H__
#define __SUPERAMCPARTICLECLUSTER_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "MCParticleList.h"
#include "SuperaBase.h"
#include "SuperaMCParticleClusterData.h"
#include "SuperaTrue2RecoVoxel3D.h"
#include "larcv/core/DataFormat/Particle.h"
#include "larcv/core/DataFormat/Voxel.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaMCParticleCluster ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaMCParticleCluster : public SuperaBase {

  public:
    /// Default constructor
    SuperaMCParticleCluster(const std::string name = "SuperaMCParticleCluster");

    /// Default destructor
    ~SuperaMCParticleCluster() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    larcv::Particle MakeParticle(const supera::LArMCParticle_t& larmcp);

    bool IsTouching(const Voxel3DMeta& meta, const VoxelSet& vs1, const VoxelSet& vs2) const;
    bool IsTouching2D(const ImageMeta& meta, const VoxelSet& vs1, const VoxelSet& vs2) const;

    //std::map<int, supera::ParticleGroup> CreateParticleGroups();
    std::map<int, supera::ParticleGroup> CreateParticleGroups();

    template <typename sed_type>
    void AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
                                 std::map<int, supera::ParticleGroup>& part_grp_v,
                                 larcv::IOManager& mgr);
    void AnalyzeSimChannel(const larcv::Voxel3DMeta& meta,
                           std::map<int, supera::ParticleGroup>& part_grp_v,
                           larcv::IOManager& mgr);
    template <typename sed_type>
    void AnalyzeFirstLastStep(const larcv::Voxel3DMeta& meta,
                              std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerTouchingLEScatter(const larcv::Voxel3DMeta& meta,
                                      std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerIonizations(std::map<int, supera::ParticleGroup>& part_grp_v);
    void ApplyEnergyThreshold(std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerConversion(std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerFamilyTouching(const larcv::Voxel3DMeta& meta,
                                   std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerTouching(const larcv::Voxel3DMeta& meta,
                             std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerTouching2D(std::map<int, supera::ParticleGroup>& part_grp_v);
    void MergeShowerDeltas(std::map<int, supera::ParticleGroup>& part_grp_v);
    void DumpHierarchy(size_t trackid,
                       const std::map<int, supera::ParticleGroup>& part_grp_v) const;
    std::vector<unsigned int> ParentTrackIDs(size_t trackid) const;
    std::vector<unsigned int> ParentShowerTrackIDs(
      size_t trackid,
      const std::map<int, supera::ParticleGroup>& part_grp_v,
      bool include_lescatter = false) const;
    void SetParticleAncestory(const larcv::Voxel3DMeta& meta3d,
                              std::map<int, supera::ParticleGroup>& part_grp_v,
                              std::vector<int>& trackid2output,
                              std::vector<int>& trackid2output2d,
                              std::set<unsigned int>& mcs_trackid_s);

  private:
    int plane_index(unsigned int cryo_id, unsigned int tpc_id, unsigned int plane_id);
    size_t SemanticPriority(size_t a, size_t b) const;
    size_t _touch_threshold, _touch_threshold2d;
    supera::MCParticleList _mcpl;
    SuperaTrue2RecoVoxel3D* _true2reco;
    std::string _ref_meta3d_cluster3d;
    std::string _ref_meta3d_tensor3d;
    std::string _ref_meta2d_tensor2d;
    std::string _output_label;
    std::string _masked_true2reco_cluster3d;
    std::string _masked_true_tensor3d;
    bool _use_sed_points;
    bool _store_dedx;
    size_t _eioni_size;
    size_t _delta_size;
    size_t _compton_size;
    //double _compton_energy;
    double _edep_threshold;
    bool _use_true_pos;
    bool _use_sed;
    bool _use_sed_lite;
    bool _use_sedlite_for_firstlast;
    bool _use_sed_for_firstlast;
    bool _check_particle_validity;
    int _projection_id;
    bool _merge_shower_delta;
    larcv::BBox3D _world_bounds;
    std::vector<std::vector<std::vector<int>>> _scan;
    std::vector<size_t> _semantic_priority;
    size_t _valid_nplanes;
    std::vector<larcv::ImageMeta> _meta2d_v;
    bool
      _useOrigTrackID; ///< Whether to use origTrackID or trackID from SimChannel/SimEnergyDeposit
  };

  /**
     \class larcv::SuperaMCParticleClusterFactory
     \brief A concrete factory class for larcv::SuperaMCParticleCluster
  */
  class SuperaMCParticleClusterProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCParticleClusterProcessFactory()
    {
      ProcessFactory::get().add_factory("SuperaMCParticleCluster", this);
    }
    /// dtor
    ~SuperaMCParticleClusterProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaMCParticleCluster(instance_name);
    }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
