# Data-ready fhicl files

These fhicl files are targeted to run Supera on data files (as opposed to MC).

Currently organized to run on both cryostats.
* `icarus_supera_data_all_cryo.fcl` is the Supera config which lives in `../include/`
* `icarus_supera_data.fcl` is the fhicl to actually run Supera.
Note that currently it will run only Supera.
It assumes the data files are at stage1 and processed appropriately already.
* `icarus_supera_data_larsoft.fcl` also runs everything else (DAQ, stage0, stage1 etc)
You should pick what you need based on your data file stage of processing!

