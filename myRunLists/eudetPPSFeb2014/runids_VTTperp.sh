#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-02-04

echo "starting script"

for GEOID in 198 199 200 201 202 203 204
do

	echo "process geoid"$GEOID
	./bdrive/bdrive_args2.sh -l myRunLists/eudetPPSFeb2014/geoid$GEOID eudetPPSFeb2014FEI4

done

echo "script finished"
