# SBND Supera Notes

## EM Shower

Keeps EM showers from being clustered with tracks

```
services.LArG4Parameters.KeepEMShowerDaughters: true # This doesn't do anything
services.ParticleListAction.keepEMShowerDaughters: false # This does since we are using legacy/refactored g4
```

Setting both to false produces not cluster_labels.
Both to true works but still produce missing ancestor and start points for tracks and photons respectively.

It turns out setting `services.LArG4Parameters.KeepEMShowerDaughters: true` is what's needed to keep cluster labels. But the other as false will merge tracks with showers.

## TPC configuration

Set TPC to [0,1] and Cryo to [0,0] when relevant (see fcls)

## Ancestor track start not set

Sedlite and sed set to true in fcl

## SuperaMCParticleCluster


Currently sets track ancestor info to dumby values (default values)

Sets photon start point to dumby (default) values

Orig track ID is -999 for all particles

After creating particle groups, all of the first_pt values are set to dumby values

Create Particle Groups makes particles if the pdg is not a nucleus, particle first step is set but not particle group in create particles function

SBND missing sed_lite time from `AnalyzeSimChannel` function. But it's not being set in there .. 

Trying to figure out if `AnalyzeFirstLastStep` fills them for ICARUS

SBND missing origTrackID - set to -999 for all parts

`sed.TrackID()` is always 0 for SBND - perhaps not filled correctly at G4 stage

`.../sbncode/LArG4` contains the `G4inforeducer` which returns 0 for all orig track IDs but is nonzero for the normal track ids.

Francois mentioned they use a different module for the `G4inforeducer` and there's an additional module in there. perhpas I'll try that if setting to the normal track ID doesn't work.

**Changing `G4inforreducer` to use trackID instead of origTrackID worked for fixing SBND's track and ancestor issue!**









