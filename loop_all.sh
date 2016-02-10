//FE-I3 50um desy small carton box
make -j4; ./setup -l myRunLists/eudetPPSNov2013/geoid100_mod -f -c eudetPPSNov2013FEI34; ./tbmon -l myRunLists/eudetPPSNov2013/geoid100_mod -f -org -c eudetPPSNov2013FEI34

//FE-I3 125 desy desy do-box
make -j4; ./setup -l myRunLists/eudetPPSNov2013/geoid165 -f -c eudetPPSNov2013FEI34;  ./tbmon -l myRunLists/eudetPPSNov2013/geoid165 -f -org -c eudetPPSNov2013FEI34

//FE-I3 50um cern do-box
make -j4; ./setup -l myRunLists/eudetPPSSep2012/geoid180 -f -c eudetPPSSep2012FEI3;  ./tbmon -l myRunLists/eudetPPSSep2012/geoid180 -f -org -c eudetPPSSep2012FEI3

//FE-I4 cern do-box


//FE-I4 desy

