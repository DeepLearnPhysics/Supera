# Job configuration examples

This folder contains many examples of job fhicl files and Supera configurations.

As a reminder, Supera runs as an *analyzer*. A Supera analyzer module needs to be
defined in your art job fhicl file. That analyzer takes several parameters. The most
important one is `supera_params` which should contain the filename of a Supera-specific
configuration. These follow a LArCV configuration format, *not a normal fhicl file format*
despite their filename ending in `.fcl`.

As a general rule of thumb,

* `run_*.fcl` are full standalone art fhicl files that can be run with `lar -c your_config.fcl`.
* `supera_*.fcl` are Supera-specific configurations written in LArCV configuration format.

The latter's filename needs to specified in the configuration of a Supera analyzer
defined in a `run_*.fcl` fhicl file.

**Note that all of the fhicl files in this folder are not necessarily maintained.
Use caution and exercise judgment. Do not blindly trust any of these configurations.**

As of 9/27/22, the following fhicl are in active usage and can be a good starting point:
* `run_mpvmpr.fcl` and `supera_mpvmpr.fcl`
* `run_nue.fcl` and `supera_nue.fcl`
* `run_corsika.fcl` and `supera_corsika.fcl`

Note that they represent a full simulation path from generator all the way to Supera.
If you are looking to run the Supera-only stage, look into the `include` folder.

## Folders
* `data` contains ready-to-run art fhicl files to run on data. Make sure to comment / uncomment
exactly which modules you need to run, depending on the processing stage of your data file.
* `test_production` should not concern anyone. They are kept for reference as they were
used to cross-check changes in larsoft when merging the ML pipeline with ICARUS production.
* `stages` has examples of each stage of the simulation pipeline broken down. If you need to
run a specific stage for debugging for example, that's a starting point.
* `include` has all the fhicl files defining configurations usable directly by ICARUS production.
This folder is meant to be more stable than anything else in this `job` folder. If you make changes
in there, please make sure to read its README.

