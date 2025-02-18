#ifndef __SUPERA_MCPARTICLEHELPER_CXX__
#define __SUPERA_MCPARTICLEHELPER_CXX__

#include "MCParticleHelper.h"
#include "FMWKInterface.h"
#include "larcv/core/Base/larbys.h"
#include <TLorentzVector.h> // ROOT
#include <set>
namespace supera {

  void MCParticleHelper::configure(const supera::Config_t& cfg)
  {
    LARCV_DEBUG() << "start" << std::endl;
    set_verbosity((::larcv::msg::Level_t)(cfg.get<unsigned short>("Verbosity", logger().level())));
    _apply_sce = cfg.get<bool>("ApplySCE");
  }

  ::larcv::Particle MCParticleHelper::MakeParticle(const supera::LArMCTrack_t& mct,
                                                   const larcv::Voxel3DMeta& meta3d) const
  {
    LARCV_DEBUG() << "start" << std::endl;
    LARCV_INFO() << "Assessing MCTrack G4Track ID = " << mct.TrackID() << " PdgCode "
                 << mct.PdgCode() << std::endl;
    double xyz[3] = {0.};

    ::larcv::Particle res;
    res.shape(::larcv::kShapeTrack);
    //res.Type(::larcv::PdgCode2ROIType(mct.PdgCode()));
    if (mct.size())
      res.energy_deposit(mct.front().E() - mct.back().E());
    else
      res.energy_deposit(0);
    res.energy_init(mct.Start().E());

    xyz[0] = mct.Start().X();
    xyz[1] = mct.Start().Y();
    xyz[2] = mct.Start().Z();

    if (_apply_sce) ApplySCE(xyz);
    res.position(xyz[0], xyz[1], xyz[2], mct.Start().T());

    xyz[0] = mct.End().X();
    xyz[1] = mct.End().Y();
    xyz[2] = mct.End().Z();
    if (_apply_sce) ApplySCE(xyz);
    res.end_position(xyz[0], xyz[1], xyz[2], mct.End().T());

    res.creation_process(mct.Process());

    if (meta3d.empty()) {

      if (mct.size() > 0) {
        auto const& first_step = mct.front();
        xyz[0] = first_step.X();
        xyz[1] = first_step.Y();
        xyz[2] = first_step.Z();
        LARCV_DEBUG() << "***--- TrackID" << mct.TrackID() << "*********------- first_step track"
                      << first_step.X() << "," << first_step.Y() << "," << first_step.Z()
                      << std::endl;
        if (_apply_sce) ApplySCE(xyz);
        res.first_step(xyz[0], xyz[1], xyz[2], first_step.T());
      }
      if (mct.size() > 1) {
        auto const& last_step = mct.back();
        LARCV_DEBUG() << "***--- TrackID" << mct.TrackID() << "*********------- last_step track"
                      << last_step.X() << "," << last_step.Y() << "," << last_step.Z() << std::endl;
        res.last_step(last_step.X(), last_step.Y(), last_step.Z(), last_step.T());
        double length = 0;
        for (size_t step_idx = 1; step_idx < mct.size(); ++step_idx) {
          auto const& step1 = mct[step_idx - 1];
          auto const& step2 = mct[step_idx];
          length += sqrt(pow(step1.X() - step2.X(), 2) + pow(step1.Y() - step2.Y(), 2) +
                         pow(step1.Z() - step2.Z(), 2));
        }
        res.distance_travel(length);
      }
    }
    else {

      int first_step = -1;
      for (size_t i = 0; i < mct.size(); ++i) {
        auto const& step = mct[i];
        auto id = meta3d.id(step.X(), step.Y(), step.Z());
        if (id == larcv::kINVALID_VOXELID) continue;
        xyz[0] = step.X();
        xyz[1] = step.Y();
        xyz[2] = step.Z();
        if (_apply_sce) ApplySCE(xyz);
        LARCV_DEBUG() << "***--- TrackID" << mct.TrackID() << "*********------- first_step track"
                      << step.X() << "," << step.Y() << "," << step.Z() << std::endl;
        res.first_step(xyz[0], xyz[1], xyz[2], step.T());
        first_step = i;
        break;
      }
      int last_step = first_step;
      for (size_t i = first_step; i < mct.size(); ++i) {
        auto const& step = mct[i];
        auto id = meta3d.id(step.X(), step.Y(), step.Z());
        if (id == larcv::kINVALID_VOXELID) break;
        xyz[0] = step.X();
        xyz[1] = step.Y();
        xyz[2] = step.Z();
        if (_apply_sce) ApplySCE(xyz);
        LARCV_DEBUG() << "***--- TrackID" << mct.TrackID() << "*********------- last_step track"
                      << step.X() << "," << step.Y() << "," << step.Z() << std::endl;
        res.last_step(xyz[0], xyz[1], xyz[2], step.T());
        last_step = i;
      }
      if (first_step > 0 && first_step != last_step) {
        double length = 0;
        for (int step_idx = first_step; step_idx < last_step; ++step_idx) {
          auto const& step1 = mct[step_idx];
          auto const& step2 = mct[step_idx + 1];
          length += sqrt(pow(step1.X() - step2.X(), 2) + pow(step1.Y() - step2.Y(), 2) +
                         pow(step1.Z() - step2.Z(), 2));
        }
        res.distance_travel(length);
      }
    }
    res.momentum(mct.Start().Px(), mct.Start().Py(), mct.Start().Pz());
    res.pdg_code(mct.PdgCode());
    res.parent_pdg_code(mct.MotherPdgCode());
    res.track_id(mct.TrackID());
    res.parent_track_id(mct.MotherTrackID());

    xyz[0] = mct.MotherStart().X();
    xyz[1] = mct.MotherStart().Y();
    xyz[2] = mct.MotherStart().Z();
    if (_apply_sce) ApplySCE(xyz);

    res.parent_position(xyz[0], xyz[1], xyz[2], mct.MotherStart().T());
    /*
    res.parent_momentum(mct.MotherStart().Px(),
                        mct.MotherStart().Py(),
                        mct.MotherStart().Pz());
    */
    LARCV_INFO() << res.dump();
    return res;
  }

  ::larcv::Particle MCParticleHelper::MakeParticle(const supera::LArMCShower_t& mcs) const
  {
    LARCV_DEBUG() << "start" << std::endl;
    LARCV_INFO() << "Assessing MCShower G4Track ID = " << mcs.TrackID() << " PdgCode "
                 << mcs.PdgCode() << std::endl;

    double xyz[3] = {0.};

    ::larcv::Particle res;
    res.shape(::larcv::kShapeShower);
    //res.Type(::larcv::PdgCode2ROIType(mcs.PdgCode()));
    res.energy_deposit(mcs.DetProfile().E());
    //res.energy_deposit(0);
    res.energy_init(mcs.Start().E());

    xyz[0] = mcs.Start().X();
    xyz[1] = mcs.Start().Y();
    xyz[2] = mcs.Start().Z();
    if (_apply_sce) ApplySCE(xyz);
    res.position(xyz[0], xyz[1], xyz[2], mcs.Start().T());

    xyz[0] = mcs.End().X();
    xyz[1] = mcs.End().Y();
    xyz[2] = mcs.End().Z();
    if (_apply_sce) ApplySCE(xyz);
    res.end_position(xyz[0], xyz[1], xyz[2], mcs.End().T());

    res.creation_process(mcs.Process());

    auto const& first_step = mcs.DetProfile();
    xyz[0] = first_step.X();
    xyz[1] = first_step.Y();
    xyz[2] = first_step.Z();
    LARCV_DEBUG() << "***--- TrackID" << mcs.TrackID() << "*********------- first_step shower"
                  << first_step.X() << "," << first_step.Y() << "," << first_step.Z() << std::endl;
    if (_apply_sce) ApplySCE(xyz);
    res.first_step(xyz[0], xyz[1], xyz[2], first_step.T());

    res.momentum(mcs.Start().Px(), mcs.Start().Py(), mcs.Start().Pz());
    res.pdg_code(mcs.PdgCode());
    res.parent_pdg_code(mcs.MotherPdgCode());
    res.track_id(mcs.TrackID());
    res.parent_track_id(mcs.MotherTrackID());

    xyz[0] = mcs.MotherStart().X();
    xyz[1] = mcs.MotherStart().Y();
    xyz[2] = mcs.MotherStart().Z();
    if (_apply_sce) ApplySCE(xyz);
    res.parent_position(xyz[0], xyz[1], xyz[2], mcs.MotherStart().T());
    /*
    res.parent_momentum(mcs.MotherStart().Px(),
                        mcs.MotherStart().Py(),
                        mcs.MotherStart().Pz());
    */
    LARCV_INFO() << res.dump();
    return res;
  }

}
#endif

// Local Variables:
// mode: c++
// End: