/**
 * \file SegLabelGhost.h
 *
 * \ingroup Supera
 *
 * \brief Class def header for a class SegLabelGhost
 *
 * @author kvtsang
 */

/** \addtogroup Supera

    @{*/
#ifndef __SEGLABELGHOST_H__
#define __SEGLABELGHOST_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
namespace larcv {

  /**
     \class ProcessBase
     User defined class SegLabelGhost ... these comments are used to generate
     doxygen documentation!
  */
  class SegLabelGhost : public ProcessBase {

  public:
    /// Default constructor
    SegLabelGhost(const std::string name = "SegLabelGhost");

    /// Default destructor
    ~SegLabelGhost() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::string _reco_prod;
    std::string _mc_prod;
    std::string _out_prod;
    size_t _min_num_voxel;
  };

  /**
     \class larcv::SegLabelGhostFactory
     \brief A concrete factory class for larcv::SegLabelGhost
  */
  class SegLabelGhostProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SegLabelGhostProcessFactory() { ProcessFactory::get().add_factory("SegLabelGhost", this); }
    /// dtor
    ~SegLabelGhostProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name)
    {
      return new SegLabelGhost(instance_name);
    }
  };

}

#endif
/** @} */ // end of doxygen group
