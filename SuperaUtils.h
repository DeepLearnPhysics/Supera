/**
 * @file SuperaUtils.h
 * @brief A collection of small classes used in various Supera modules.
 * @author Patrick Tsang
 * @version v0.1.0
 * @date 2021-09-01
 */

#ifndef __SUPERAUTILS__
#define __SUPERAUTILS__

#include "larcv/core/DataFormat/Voxel.h"
#include <set>

namespace supera {
  namespace utils {

    /**
     * @brief A faster implementation of larcv::VoxelSet.
     *
     * What's wrong with larcv::VoxelSet?
     * ------------------------------------------
     * The underlying container is in std::vector.
     * O(N) for random insert to i-th element
     *   - extend container size by 1
     *   - shift elements by 1 start from i
     *   - insert new element to location i
     * FastVoxelSet use std::set for storage,
     * which is O(ln(N)) for insertion.
     *
     * Why not modifty larcv::VoxelSet?
     * ------------------------------------------
     * Ideally we should. However it may break
     * analyzer's code for reading out original
     * larcv data. So FactVoxelSet is a proxy at
     * production stage to speedup insertion time
     * of a large number of voxels. The final
     * product is coverted to larcv::VoxelSet.
     */
    class FastVoxelSet {
    public:
      FastVoxelSet() { ; }

      /**
         * @brief Insert a voxel to the set
         *
         * @param id Voxel id, usually produced by larcv::Meta
         * @param value value of the voxel
         * @param add if true, add value to the exisiting voxel
         */
      void emplace(larcv::VoxelID_t id, float value, bool add)
      {
        larcv::Voxel v(id, value);
        auto itr = _voxel_set.find(v);
        if (add && itr != _voxel_set.end()) {
          v += itr->value();
          _voxel_set.erase(itr);
        }
        _voxel_set.insert(std::move(v));
      }

      /**
         * @brief Move all contents to a larcv::VoxelSet
         *
         * @param v_set Target larcv::VoxelSet
         */
      void move_to(larcv::VoxelSet& v_set)
      {
        size_t n = _voxel_set.size();
        v_set.reserve(n);
        for (auto& v : _voxel_set) {
          v_set.insert(v);
        }
        _voxel_set.clear();
      }

      /**
         * @brief size of voxel set
         *
         * @return size_t
         */
      size_t size() const { return _voxel_set.size(); }

    private:
      std::set<larcv::Voxel> _voxel_set;
    };
  }
}
#endif
