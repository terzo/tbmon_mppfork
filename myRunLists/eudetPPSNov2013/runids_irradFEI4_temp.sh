#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-02-04

echo "starting script"

for GEOID in 153 154 155 156 157 158 159 160 161 162
do

	echo "process geoid"${GEOID}
	./bdrive_temp/bdrive_args2.sh -l myRunLists/eudetPPSNov2013/geoid${GEOID} eudetPPSNov2013FEI4

done


echo "script finished"
