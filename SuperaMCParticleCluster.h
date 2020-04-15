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
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "MCParticleList.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
#include "larcv/core/DataFormat/Voxel.h"
#include "larcv/core/DataFormat/Particle.h"
#include "SuperaMCParticleClusterData.h"
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

    std::vector<supera::ParticleGroup> CreateParticleGroups();

    void AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
				 std::vector<supera::ParticleGroup>& part_grp_v,
         larcv::IOManager& mgr);
    void AnalyzeSimChannel(const larcv::Voxel3DMeta& meta,
			   std::vector<supera::ParticleGroup>& part_grp_v,
			   larcv::IOManager& mgr);
    void AnalyzeFirstLastStep(const larcv::Voxel3DMeta& meta,
			      std::vector<supera::ParticleGroup>& part_grp_v);
    void MergeShowerIonizations(std::vector<supera::ParticleGroup>& part_grp_v);
    void ApplyEnergyThreshold(std::vector<supera::ParticleGroup>& part_grp_v);
    void MergeShowerConversion(std::vector<supera::ParticleGroup>& part_grp_v);
    void MergeShowerFamilyTouching(const larcv::Voxel3DMeta& meta,
				   std::vector<supera::ParticleGroup>& part_grp_v);
    void MergeShowerTouching(const larcv::Voxel3DMeta& meta,
			     std::vector<supera::ParticleGroup>& part_grp_v);
    void MergeShowerTouching2D(std::vector<supera::ParticleGroup>& part_grp_v);

    void MergeShowerDeltas(std::vector<supera::ParticleGroup>& part_grp_v);
    void DumpHierarchy(size_t trackid,
		       const std::vector<supera::ParticleGroup>& part_grp_v) const;
    std::vector<unsigned int> ParentTrackIDs(size_t trackid) const;
  private:
    int plane_index(unsigned int cryo_id, unsigned int tpc_id, unsigned int plane_id);
    size_t SemanticPriority(size_t a, size_t b) const;
    supera::MCParticleList _mcpl;
    std::string _ref_meta3d_cluster3d;
    std::string _ref_meta3d_tensor3d;
    std::string _ref_meta2d_tensor2d;
    std::string _output_label;
    std::string _masked_true2reco_cluster3d;
    std::string _masked_true_tensor3d;
    bool _use_sed_points;
    size_t _eioni_size;
    size_t _delta_size;
    size_t _compton_size;
    double _compton_energy;
    double _edep_threshold;
    bool _use_true_pos;
    bool _use_sed;
    bool _check_particle_validity;
    int  _projection_id;
    larcv::BBox3D _world_bounds;
    std::vector<std::vector<std::vector<int> > > _scan;
    std::vector<size_t> _semantic_priority;
    size_t _valid_nplanes;
    std::vector<larcv::ImageMeta> _meta2d_v;
  };

  /**
     \class larcv::SuperaMCParticleClusterFactory
     \brief A concrete factory class for larcv::SuperaMCParticleCluster
  */
  class SuperaMCParticleClusterProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCParticleClusterProcessFactory() { ProcessFactory::get().add_factory("SuperaMCParticleCluster", this); }
    /// dtor
    ~SuperaMCParticleClusterProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaMCParticleCluster(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
