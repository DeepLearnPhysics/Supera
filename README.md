# Supera
Supera is designed to slice up the data products output from larsoft and package them to be consumed by the mlreco workflow. The bulk of the work is done by `SuperaMCParticleCluster` which consumes inputs from `simenergydeposit(lite)`, cluster3D, and MC truth products such as `MCTrack`, `MCShower`, and `MCParticle`.

## Verbosity
Verbosity can be set for any module used in supera by setting `Verbosity : n`

* `n = 0` shows `LARCV_DEBUG` and above messages
* `n = 1` shows `LARCV_INFO` and above messages
* `n = 2` shows `LARCV_NORMAL` and above messages
* `n = 3` shows `LARCV_WARNING` and above messages
* `n = 4` shows `LARCV_ERROR` and above messages
* `n = 5` only shows `LARCV_CRITICAL` and messages

## SuperaMCParticleCluster
This class is responsible for managing Monte Carlo (MC) particle clusters within the LArCV framework. The class is authored by drinkkazu.

### Features
* It provides functions to configure, initialize, process, and finalize instances of the class.
* It allows creation of particles and particle groups.
* It includes methods for analyzing simulated energy deposits and channels.
* It has functions for merging shower touchings, conversions, ionizations, and deltas.
* It includes functions for merging shower family touchings.
* The class also contains private variables and methods to work with particle clusters.
### Dependencies
The class depends on several other modules and classes in the LArCV framework, including:

* `SuperaBase`
* `SuperaTrue2RecoVoxel3D`
* `FMWKInterface`
* `MCParticleList`
* Various data formats from `larcv/core/DataFormat`
* `SuperaMCParticleClusterData`

The class is built upon a factory pattern, with the corresponding `SuperaMCParticleClusterProcessFactory` class provided for creating instances of `SuperaMCParticleCluster`.

# Building Supera for ProtoDUNE

In general, **Supear** works fine with latest larsoft release on _cvmfs_ (linux only, sorry for mac users).
The following has been tested with the tags specified. It depends on [[larcv2][https://github.com/DeepLearnPhysics/larcv2]].

**For building**
```
# ===========================
# Prepare a working directory
# ===========================
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup larsoft v08_21_00 -q e17:prof

mkdir MyWorkDir && cd MyWorkDir
mrb newDev
source MyWorkDir/localProducts_larsoft_v08_21_00_e17_prof/setup

# ===============
# checkout larcv2
# ===============
cd $MRB_SOURCE
git clone https://github.com/DeepLearnPhysics/larcv2
cd larcv2
source configure.sh
make

# ================
# checkout dunetpc
# ================
cd $MRB_SOURCE
mrb g dunetpc

# ===============
# checkout supera
# ===============
cd $MRB_SOURCE/dunetpc/dune
git clone -b pdune https://github.com/kvtsang/Supera.git
cd Supera
source setup.sh pdune

# ============================================
# edit $MRB_SOURCE/dunetpc/dune/CMakeLists.txt
# add this line "add_subdirectory(Supera)"
# ============================================

cd $MRB_BUILDDIR
mrbsetenv
mrb i             # or mrb i -j <n> for parallel build
```

**For running**
```
# ==========================================
# Do it once, before executing "lar" command
# ==========================================
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source MyWorkDir/localProducts_larsoft_v08_21_00_e17_prof/setup
source MyWorkDir/larcv2/configure.sh
```
