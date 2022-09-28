# Data-ready fhicl files

These fhicl files are targeted to run Supera on data files (as opposed to MC).

Currently organized to run on each cryostat separately:
* `supera_data_cryoE.fcl` and `supera_data_cryoW.fcl` for 2 separate analyzer modules
(one per cryostat), the configs lives in `../include/`
* `icarus_supera_data.fcl` is the fhicl to actually run Supera.
Note that currently it will run only Supera but everything else (DAQ, stage0, stage1 etc)
is just commented out, i.e. it
assumes the data files are at stage1 and processed appropriately already.
You should uncomment what you need based on your data file stage of processing!

**Important note**
For CryoW, the points are shifted along X-axis to overlap with CryoE layout by SuperaSpacePoint module.

The output will be 2 LArCV files (`larcv_cryoE.root` and `larcv_cryoW.root`).
