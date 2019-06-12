# Building Supera for ProtoDUNE

In general, **Supear** works fine with latest larsoft release on _cvmfs_ (linux only, sorry for mac users).
The following has been tested with the tags specified. It depends on [[larcv2][https://github.com/DeepLearnPhysics/larcv2]].

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
