#! /usr/local/bin/bash
# run analysis of all november single chip data from november
# including the 4 calibration steps (hotpixel, align, etacorr, align2)

# VTT-NP1-4-E5 unirradiated FE-I3
./myRunLists/eudetPPSNov2013/runids_unirradFEI3.sh | tee runall_full.log

# VTT-NP1-4-E3 5e15 FE-I3
./myRunLists/eudetPPSNov2013/runids_irradFEI3.sh | tee -a runall_full.log

# VTT-NP1-8-E4 5e15 FE-I4
./myRunLists/eudetPPSNov2013/runids_irradFEI4_updated.sh | -a runall_full.log

