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

Module sets ancestor and parent info for particles

Currently sets track ancestor info to dumby values (default values)

Sets photon start point to dumby (default) values

Orig track ID is -999 for all particles



