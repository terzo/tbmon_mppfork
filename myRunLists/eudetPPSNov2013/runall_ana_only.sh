#! /usr/local/bin/bash
# run analysis of all november single chip data from november
# including the 4 calibration steps (hotpixel, align, etacorr, align2)

cp bdrive/driver_original.cc driver.cc
make

# VTT-NP1-4-E5 unirradiated FE-I3
./tbmon -l myRunLists/eudetPPSNov2013/geoid100 -f -org -c eudetPPSNov2013FEI34 | tee runall_ana_only.log

# VTT-NP1-4-E3 5e15 FE-I3
for GEOID in 163 164 165 166 167
do
	./tbmon -l myRunLists/eudetPPSNov2013/geoid${GEOID} -f -org -c eudetPPSNov2013FEI34 | tee runall_ana_only.log
done

# VTT-NP1-8-E4 5e15 FE-I4
for GEOID in 153 154 155 156 157 158 159 160 161 162
do
	./tbmon -l myRunLists/eudetPPSNov2013/geoid${GEOID} -f -org -c eudetPPSNov2013FEI4 | tee runall_ana_only.log
done
