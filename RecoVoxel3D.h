/**
 * \file RecoVoxel3D.h
 *
 * \ingroup Supera
 * 
 * \brief Class def header for a class RecoVoxel3D
 *
 * @author kvtsang
 */

/** \addtogroup Supera

    @{*/
#ifndef SUPERA_RECOVOXEL3D_H
#define SUPERA_RECOVOXEL3D_H

#include <functional>
#include "larcv/core/DataFormat/DataFormatTypes.h"

namespace larcv{
  /**
     \class RecoVoxel3D
     User defined class RecoVoxel3D... these comments are used to generate
     doxygen documentation!
  */
  class RecoVoxel3D{
    
  public:
    
    /// Default constructor
    RecoVoxel3D(VoxelID_t id, double charge=0) : _voxel_id(id), _charge(charge) {}
    
    /// Default(destructor
    ~RecoVoxel3D(){}

    VoxelID_t get_id() const { return _voxel_id; } 
    double get_charge() const { return _charge; } 

    void set_charge(double charge, bool is_add=false) {
      if (is_add)
        _charge += charge;
      else
        _charge = charge;
    }
    
    friend inline bool operator< (const RecoVoxel3D& lhs, const RecoVoxel3D& rhs) { 
      return lhs._voxel_id < rhs._voxel_id;
    }

    friend inline bool operator== (const RecoVoxel3D& lhs, const RecoVoxel3D& rhs) { 
      return lhs._voxel_id == rhs._voxel_id;
    }


  private:
    VoxelID_t _voxel_id; // voxel id
    double _charge;      // reco charge
  };
}

namespace std {
  template<>
  class hash<larcv::RecoVoxel3D> {
    public:
      size_t operator()(const larcv::RecoVoxel3D& pt) const {
        return std::hash<larcv::VoxelID_t>()(pt.get_id());
      }
  };
}


#endif
/** @} */ // end of doxygen group 

