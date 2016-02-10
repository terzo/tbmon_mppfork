#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-02-04

echo "starting script"

for GEOID in 219 220 221 222 223
do

	echo "process geoid"$GEOID
	./bdrive/bdrive_args2.sh -l myRunLists/eudetPPSFeb2014/geoid$GEOID eudetPPSFeb2014FEI4

done

echo "script finished"
