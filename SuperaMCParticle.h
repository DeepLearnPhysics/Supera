/**
 * \file SuperaMCParticle.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaMCParticle
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMCPARTICLE_H__
#define __SUPERAMCPARTICLE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "MCParticleHelper.h"
#include "MCParticleTree.h"
#include "SuperaBase.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaMCParticle ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaMCParticle : public SuperaBase {

  public:
    /// Default constructor
    SuperaMCParticle(const std::string name = "SuperaMCParticle");

    /// Default destructor
    ~SuperaMCParticle() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    //bool StoreG4SecondaryParticle() const { return _store_g4_secondary_roi; }

    //bool StoreG4PrimaryParticle() const { return _store_g4_primary_roi; }

    //const std::vector<std::pair<size_t, size_t> >& Particle2MCNode(int roi_index) const;

    //const std::vector<std::vector<std::pair<size_t, size_t> > >&
    //Particle2MCNode() const { return _roi2mcnode_vv; }

    const supera::MCParticleTree& ParticleTree() const { return _mcpt; }

  private:
    //bool _store_part;
    //bool _store_g4_primary_part;
    //bool _store_g4_secondary_part;
    //bool _store_pixel2d;
    //bool _store_voxel3d;
    std::vector<larcv::Particle> _part_v;
    supera::MCParticleTree _mcpt;
    supera::MCParticleHelper _mcpart_maker;
    //std::vector<std::vector<std::pair<size_t, size_t> > > _part2mcnode_vv;

    unsigned short _pass_origin;
    std::vector<int> _filter_pdg;
    std::vector<double> _filter_min_einit;
    std::vector<double> _filter_min_edep;

    double _shower_min_einit;
    double _shower_min_edep;

    double _track_min_einit;
    double _track_min_edep;

    std::string _output_label;
    std::string _ref_meta3d_cluster3d;
    std::string _ref_meta3d_tensor3d;

    //size_t _filter_min_cols;
    //size_t _filter_min_rows;

    bool FilterNode(const supera::MCNode& node) const;
    larcv::Particle MakeParticle(const supera::MCNode& node, const larcv::Voxel3DMeta& meta) const;
  };

  /**
     \class larcv::SuperaMCParticleFactory
     \brief A concrete factory class for larcv::SuperaMCParticle
  */
  class SuperaMCParticleProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMCParticleProcessFactory()
    {
      ProcessFactory::get().add_factory("SuperaMCParticle", this);
    }
    /// dtor
    ~SuperaMCParticleProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SuperaMCParticle(instance_name);
    }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group
