#ifndef __FMWKINTERFACE_H__
#define __FMWKINTERFACE_H__

//#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcore/Geometry/Geometry.h"
//#include "fhiclcpp/ParameterSet.h"
#include "larcv/core/Base/PSet.h"
#include "larcv/core/DataFormat/DataFormatTypes.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "lardataobj/MCBase/MCMiniPart.h"
#include "ExperimentTypes.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/MCBase/MCParticleLite.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimEnergyDepositLite.h"

namespace supera {
  typedef larcv::PSet Config_t;
  typedef recob::Wire LArWire_t;
  typedef raw::OpDetWaveform LArOpDigit_t;
  typedef recob::Hit LArHit_t;
  typedef simb::MCTruth LArMCTruth_t;
  typedef simb::MCParticle LArMCParticle_t;
  typedef sim::MCParticleLite LArMCMiniPart_t;
  typedef sim::MCTrack LArMCTrack_t;
  typedef sim::MCShower LArMCShower_t;
  typedef sim::SimChannel LArSimCh_t;
  typedef sim::MCStep LArMCStep_t;
  typedef sim::SimEnergyDeposit LArSimEnergyDeposit_t;
  typedef sim::SimEnergyDepositLite LArSimEnergyDepositLite_t;
  typedef recob::SpacePoint LArSpacePoint_t;
  typedef recob::OpFlash LArOpFlash_t;
}
//
// Utility functions (geometry, lar properties, etc.)
//
namespace supera {

  //typedef ::fhicl::ParameterSet Config_t;

  //
  // LArProperties
  //

  /// DriftVelocity in cm/us
  double DriftVelocity();

  //
  // Geometry
  //

  /// Channel number to wire ID
  ::geo::WireID ChannelToWireID(unsigned int ch);

  /// Channel number to image column number
  double ChannelToImageX(unsigned int ch);

  /// Channel number to projection ID
  larcv::ProjectionID_t ChannelToProjectionID(unsigned int ch);

  /// Number of channels
  unsigned int Nchannels();
  /*
  /// Number of planes (deprecated)
  unsigned int Nplanes();

  /// Number of projections
  unsigned int NProjections();

  /// Number of wires
  unsigned int Nwires(unsigned int plane);

  /// Nearest wire
  unsigned int NearestWire(const TVector3& xyz, unsigned int plane);

  /// Nearest wire
  unsigned int NearestWire(const double* xyz, unsigned int plane);

  /// Angle from z-axis
  double WireAngleToVertical(unsigned int plane);

  /// Wire pitch
  double WirePitch(size_t plane);

  /// Detector height
  double DetHalfHeight();

  /// Detector width
  double DetHalfWidth();

  /// Detector length
  double DetLength();
  */
  //
  // DetectorClockService
  //

  /// Number of time ticks
  unsigned int NumberTimeSamples();

  /// G4 time to TPC tick
  int TPCG4Time2Tick(double ns);

  /// G4 time to TPC tick
  int TPCG4Time2TDC(double ns);

  /// per-plane tick offset
  double PlaneTickOffset(size_t plane0, size_t plane1);

  /// TPC TDC to Tick
  double TPCTDC2Tick(double tdc);

  /// Trigger time [us] offset
  double TriggerOffsetTPC();

  /// TPC sampling period
  double TPCTickPeriod();

  /// Trigger time
  double TriggerTime();

  //
  // SpaceChargeService
  //
  /// Truth position to shifted
  void ApplySCE(double& x, double& y, double& z);
  /// Truth position to shifted
  void ApplySCE(double* xyz);
  /// Truth position to shifted
  void ApplySCE(TVector3& xyz);

}

#endif
