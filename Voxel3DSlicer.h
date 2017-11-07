/**
 * \file Voxel3DSlicer.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class Voxel3DSlicer
 *
 * @author kazuhiro
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __VOXEL3DSLICER_H__
#define __VOXEL3DSLICER_H__
#include "FMWKInterface.h"
#include "WireRange3D.h"
#include "ImageMetaMakerBase.h"

namespace supera {

  /**
     \class Voxel3DSlicer
     User defined class Voxel3DSlicer ... these comments are used to generate
     doxygen documentation!
  */
  class Voxel3DSlicer : public ImageMetaMakerBase {

  public:

    /// Default constructor
    Voxel3DSlicer() : ImageMetaMakerBase("Voxel3DSlicer")
      , _slicer()
    {}

    /// Default destructor
    ~Voxel3DSlicer() {}

    static bool Is(const ImageMetaMakerBase* ptr)
    { return (ptr && ptr->name() == "Voxel3DSlicer"); }

    void configure(const supera::Config_t&);

    void AddConstraint(double x, double y, double z);

    void AddConstraint(const supera::LArMCTruth_t& mctruth);

    void AddConstraint(const std::vector<supera::LArMCTruth_t>& mctruth_v);

    void
    GenerateMeta(const std::vector<supera::LArSimCh_t>& simch_v,
                 const int time_offset);

    void
    GenerateMeta(const std::vector<supera::LArSimCh_t>& simch_v,
                 const int time_offset,
                 const std::vector<int>& trackid_v);

    void
    GenerateMeta(const int time_offset);

    void ClearEventData();

    bool Test() const;

    inline const std::vector<larcv::ImageMeta>& Meta() const
    { return _meta_v;}

    inline const larcv::Voxel3DMeta& Meta3D() const
    { return _meta3d; }

  private:

    unsigned short _origin;
    supera::WireRange3D _slicer;
    bool _apply_sce;
    double _t0_g4ns;
    std::vector<larcv::ImageMeta> _meta_v;
    larcv::Voxel3DMeta _meta3d;

    void DeriveMeta(std::vector<larcv::ImageMeta>& meta_v,
		    larcv::Voxel3DMeta& meta3d,
		    const std::vector<supera::GridPoint3D>& point_v,
		    const int time_offset ) const;

    void
    GenerateMeta(const std::vector<supera::LArSimCh_t>& simch_v,
                 const int time_offset,
                 const std::vector<int>& trackid_v,
                 const bool use_track_id);
  };

}
#endif
/** @} */ // end of doxygen group

