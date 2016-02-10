#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-04-01

echo "starting script"

for GEOID in 205 206 207 208 209 210
do

	echo "process geoid"$GEOID
	./bdrive/bdrive_args2.sh -l myRunLists/eudetPPSFeb2014/geoid$GEOID eudetPPSFeb2014FEI4

done

echo "script finished"
