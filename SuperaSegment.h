/**
 * \file SuperaSegment.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class SuperaSegment
 *
 * @author josh
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERASEGMENT_H__
#define __SUPERASEGMENT_H__
#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "ImageMetaMaker.h"
#include "ParamsImage2D.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include <string>

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaSegment ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaSegment : public SuperaBase,
    public supera::ParamsImage2D,
    public supera::ImageMetaMaker {

  public:

    /// Default constructor
    SuperaSegment(const std::string name="SuperaSegment");

    /// Default destructor
    ~SuperaSegment(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    double getTick( const std::vector<double>& step, const double trig_time=4050.0);

    double getTick( const supera::LArMCStep_t& step, const double trig_time=4050.0);

    std::vector<int> getFirstStepPosInsideImage( const supera::LArMCTrack_t& track, const larcv::ImageMeta& meta, const double trig_time,
						 const bool startAtstart, const double max_step_size, const double fv_border); //, const double fv_border  another argument to function

    std::vector<int> getProjectedImagePixel( const std::vector<double>& pos3d, const larcv::ImageMeta& meta, const int nplanes, const double fracpixborder = 1.5 );

    int doesTrackCrossImageBoundary( const supera::LArMCTrack_t& track, const larcv::ImageMeta& meta, const double trig_time);



    float dwall( const std::vector<float>& pos, int& boundary_type );

    double dwall( const std::vector<double>& pos, int& boundary_type );


    //std::vector<int> getImageBoundaryCrossingPoint( const supera::LArMCTrack_t& track, std::vector<double>& crossingpt, const larcv::ImageMeta& meta,
    //  const double boundary_tick_buffer, const double trig_time);

    //int doesTrackCrossImageBoundary( const supera::LArMCTrack_t& track, const larcv::ImageMeta& meta, const double trig_time);

    // double getTick( const std::vector<double>& step, const double trig_time=4050.0);
    //
    // double getTick( const supera::LArMCStep_t& step, const double trig_time=4050.0);
    //
    //std::vector<int> getProjectedImagePixel( const std::vector<double>& pos3d, const larcv::ImageMeta& meta, const int nplanes, const double fracpixborder = 1.5 );

    unsigned short _origin;
    std::string m_ancestor_label;
    std::string m_instance_label;
    std::string m_wire_label;
  };


  /**
     \class larcv::SuperaSegmentFactory
     \brief A concrete factory class for larcv::SuperaSegment
  */
  class SuperaSegmentProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaSegmentProcessFactory() { ProcessFactory::get().add_factory("SuperaSegment",this); }
    /// dtor
    ~SuperaSegmentProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string segment_name) { return new SuperaSegment(segment_name); }
  };


  }
#endif
/** @} */ // end of doxygen group
