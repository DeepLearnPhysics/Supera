/**
 * @file SpacePointHits.h
 * @brief Hits associated with each space point, without making art copy the association.
 */

#ifndef __SUPERA_SPACEPOINTHITS_H__
#define __SUPERA_SPACEPOINTHITS_H__

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include <vector>

namespace supera {
  namespace utils {

    /**
     * @brief The hits of each space point, in the order they appear in the association.
     *
     * This is what art::FindManyP<recob::Hit>(space points, ...) returns. Cluster3D stores
     * art::Assns<recob::Hit, recob::SpacePoint>; FindManyP asks for the opposite orientation, which
     * makes art build a full copy of the association (120 MB on a busy SBND event) that then stays
     * in the event. Reading the association as stored avoids the copy.
     */
    class SpacePointHits {
    public:
      SpacePointHits(art::ValidHandle<std::vector<recob::SpacePoint>> const& space_points,
                     art::Event const& event,
                     art::InputTag const& tag)
        : _hits(space_points->size())
      {
        auto const& assns =
          *event.getValidHandle<art::Assns<recob::Hit, recob::SpacePoint>>(tag);
        auto const id = space_points.id();
        for (auto const& hit_sp : assns) {
          auto const& sp = hit_sp.second;
          if (sp.id() == id && sp.key() < _hits.size()) _hits[sp.key()].push_back(hit_sp.first);
        }
      }

      /// Same use as art::FindManyP::get() on an empty vector
      void get(size_t i, std::vector<art::Ptr<recob::Hit>>& hits) const { hits = _hits.at(i); }

    private:
      std::vector<std::vector<art::Ptr<recob::Hit>>> _hits;
    };
  }
}

#endif
