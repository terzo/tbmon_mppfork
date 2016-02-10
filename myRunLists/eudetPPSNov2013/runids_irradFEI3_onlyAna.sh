#! /bin/bash
# script needs to be called from trunk directory

echo "starting script"
echo ">>> execute make"
make

echo "process geoid163"
echo ">>> execute tbmon"
./tbmon -l myRunLists/eudetPPSNov2013/geoid163 -f -org -c eudetPPSNov2013FEI34

echo "process geoid164"
echo ">>> execute tbmon"
./tbmon -l myRunLists/eudetPPSNov2013/geoid164 -f -org -c eudetPPSNov2013FEI34

echo "process geoid165"
echo ">>> execute tbmon"
./tbmon -l myRunLists/eudetPPSNov2013/geoid165 -f -org -c eudetPPSNov2013FEI34

echo "process geoid166"
echo ">>> execute tbmon"
./tbmon -l myRunLists/eudetPPSNov2013/geoid166 -f -org -c eudetPPSNov2013FEI34

echo "process geoid167"
echo ">>> execute tbmon"
./tbmon -l myRunLists/eudetPPSNov2013/geoid167 -f -org -c eudetPPSNov2013FEI34 

echo "script finished"
