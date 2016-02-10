#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-03-19

echo "starting script"
echo ">>> execute make"
make

for GEOID in 153 154 155 156 157 158 159 160 161 162
do

	echo "process geoid"${GEOID}
	echo ">>> execute tbmon"
	./tbmon -l myRunLists/eudetPPSNov2013/geoid${GEOID}_updated eudetPPSNov2013FEI4

done


echo "script finished"