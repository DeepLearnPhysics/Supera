# FHICL files for production test

This contains fhicl files used to test the MC production workflow.
It should be run in order:
* genie
* g4
* detsim
* mcreco
* stage0
* stage1

It is based on a Summer 2022 report of what fhicl files are typically used
by the ICARUS production team for each stage, according to Francois.

In particular, it was used with `icaruscode` up to `v09_58_00` and
to test the changes made to lardataobj/larsim to merge the ML pipeline
with ICARUS official production.

Note that it might become obsolete moving forward if ICARUS production changes its
organization. Feel free to update.
