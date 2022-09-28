# Fhicl files to be used by production

This folder defines some standard fhicl files that can be used by ICARUS production.

By default it is set to use `sim::SimEnergyDepositLite` and `sim::MCParticleLite`
in order to be compatible with fhicl parameters such as `KeepEMShowerDaughters: false`
in LArG4 or the drop of `sim::SimEnergyDeposit` objects.

## Note for any developer/user

**If you make changes to one of the fhicl, please propagate to all of them as needed.**

The main fhicl to be shared/included in `icaruscode` or `sbndcode` is `supera_modules.fcl`
which defines all the analyzer modules that are made available in this folder. Be aware
that which analyzer you choose to use may depend on what your generator is (`generator` only,
`cosmgen` only or both at the same time) and MC vs data.

## Note for Supera developers

_Note that this is not ideal - lots of repetition to accomodate different setups.
As far as I can tell LArCV configuration does not let you include recursively and
modify configurations. Hence the brazen copy-paste. If you think of a better solution
or if you improve larcv2 to allow configuration includes, please update this folder._
