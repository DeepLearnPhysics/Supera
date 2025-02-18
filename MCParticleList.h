/**
 * \file MCParticleList.h
 *
 * \ingroup Supera
 *
 * \brief Class def header for a class MCParticleList
 *
 * @author kterao
 */

/** \addtogroup Supera

    @{*/
#ifndef SUPERA_MCPARTICLELIST_H
#define SUPERA_MCPARTICLELIST_H

#include "nusimdata/SimulationBase/MCParticle.h"
#include <iostream>

namespace supera {

  /**
     \class MCParticleList
     User defined class MCParticleList ... these comments are used to generate
     doxygen documentation!
  */
  class MCParticleList {

  public:
    /// Default constructor
    MCParticleList() {}

    /// Default destructor
    ~MCParticleList() {}

    void Update(const std::vector<simb::MCParticle>& larmcp_v, int run, int event);

    const std::vector<int>& PdgCode() const { return _pdgcode_v; }
    const std::vector<int>& ParentIndex() const { return _parent_index_v; }
    const std::vector<int>& ParentTrackId() const { return _parent_trackid_v; }
    const std::vector<int>& ParentPdgCode() const { return _parent_pdg_v; }
    const std::vector<int>& TrackIdToIndex() const { return _trackid2index; }
    const std::vector<int>& AncestorIndex() const { return _ancestor_index_v; }
    const std::vector<int>& AncestorTrackId() const { return _ancestor_trackid_v; }

  private:
    int _run;                         ///< a unique event identifier (together with _event)
    int _event;                       ///< a unique event identifier (together with _run)
    std::vector<int> _trackid_v;      ///< Track ID, index = std::vector<simb::MCParticle> index
    std::vector<int> _pdgcode_v;      ///< PDG code, index = std::vector<simb::MCParticle> index
    std::vector<int> _parent_index_v; ///< Parent index, index = std::vector<simb::MCParticle> index
    std::vector<int>
      _parent_trackid_v; ///< Parent track ID, index = std::vector<simb::MCParticle> index
    std::vector<int> _parent_pdg_v; ///< Parent PDG, index = std::vector<simb::MCParticle> index
    std::vector<int>
      _ancestor_index_v; ///< Ancestor index, index = std::vector<simb::MCParticle> index
    std::vector<int>
      _ancestor_trackid_v; ///< Ancestor track ID, index = std::vector<simb::MCParticle> index
    std::vector<int> _trackid2index; ///< TrackID => std::vector<simb::MCParticle> index converter
  };
}

#endif
/** @} */ // end of doxygen group
