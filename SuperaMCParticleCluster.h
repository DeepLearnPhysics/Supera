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
    void FillParticleGroups(const std::vector<simb::MCParticle>& larmcp_v,
			    std::vector<supera::ParticleGroup>& shower_grp_v,
			    std::vector<supera::ParticleGroup>& track_grp_v);
      
    void AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
				 std::vector<supera::ParticleGroup>& shower_grp_v,
				 std::vector<supera::ParticleGroup>& track_grp_v);
    void AnalyzeSimChannel(const larcv::Voxel3DMeta& meta,
			   std::vector<supera::ParticleGroup>& shower_grp_v,
			   std::vector<supera::ParticleGroup>& track_grp_v);
    void MergeShowerIonizations(std::vector<supera::ParticleGroup>& shower_grp_v);
    void ApplyEnergyThreshold(std::vector<supera::ParticleGroup>& shower_grp_v,
			      std::vector<supera::ParticleGroup>& track_grp_v);
    void MergeShowerConversion(std::vector<supera::ParticleGroup>& shower_grp_v);
    void MergeShowerTouching(const larcv::Voxel3DMeta& meta,
			     std::vector<supera::ParticleGroup>& shower_grp_v);
    void MergeShowerDeltas(std::vector<supera::ParticleGroup>& shower_grp_v,
			   std::vector<supera::ParticleGroup>& track_grp_v);
  private:
    supera::MCParticleList _mcpl;
    std::string _ref_meta3d_cluster3d;
    std::string _ref_meta3d_tensor3d;
    std::string _output_label;
    size_t _eioni_size;
    size_t _delta_size;
    size_t _compton_size;
    double _compton_energy;
    double _edep_threshold;
    bool _use_true_pos;
    bool _use_sed;
    larcv::ProjectionID_t _projection_id;
    larcv::BBox3D _world_bounds;
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

