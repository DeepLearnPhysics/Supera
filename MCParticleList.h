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
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace supera {

  /**
     \class TrackIdIndex
     Track ID => index lookup. It reads like the std::vector<int> it replaces (size() and operator[],
     with -1 for a track ID that is not there), but only stores the track IDs that exist: track IDs can
     carry large offsets, and a vector indexed by track ID was 160 MB in SBND.
  */
  class TrackIdIndex {
  public:
    /// One past the largest track ID seen so far (never shrinks, like the vector it replaces)
    size_t size() const { return _size; }
    int operator[](size_t trackid) const
    {
      auto it = _index.find(static_cast<int>(trackid));
      return it == _index.end() ? -1 : it->second;
    }
    /// Checked access: throws for a track ID at or beyond size(), like std::vector::at
    int at(size_t trackid) const
    {
      if (trackid >= _size) throw std::out_of_range("TrackIdIndex::at");
      return (*this)[trackid];
    }
    /// Forgets all entries and makes sure size() is at least min_size
    void reset(size_t min_size)
    {
      _index.clear();
      if (_size < min_size) _size = min_size;
    }
    void set(int trackid, int index)
    {
      _index[trackid] = index;
      if (trackid >= 0 && static_cast<size_t>(trackid) + 1 > _size) _size = static_cast<size_t>(trackid) + 1;
    }

  private:
    std::unordered_map<int, int> _index;
    size_t _size = 0;
  };

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
    const TrackIdIndex& TrackIdToIndex() const { return _trackid2index; }
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
    TrackIdIndex _trackid2index; ///< TrackID => std::vector<simb::MCParticle> index converter
  };
}

#endif
/** @} */ // end of doxygen group
