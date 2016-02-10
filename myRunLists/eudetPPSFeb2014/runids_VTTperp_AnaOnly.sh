#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-04-01

# modified 2014-05-13 to set MaskingRadius to 1

# file for logging output of this script
LOGFILE=myRunLists/eudetPPSFeb2014/runids_CIS2_W07_NOMOD1_AnaOnly.log



echo
echo "#############################################################"
echo "starting script" | tee $LOGFILE
echo
echo "output will be logged in "$LOGFILE
echo
echo "copy driver_original_test2.cc from bdrive folder"
echo
cp bdrive/driver_original_test2.cc driver.cc
echo "call make"
make

for GEOID in 198 199 200 201 202 203 204
do

	echo "###################################################################" \
		| tee -a $LOGFILE
	echo "process geoid"$GEOID | tee -a $LOGFILE
	echo "###################################################################" \
		| tee -a $LOGFILE
	
	./tbmon -l myRunLists/eudetPPSFeb2014/geoid$GEOID \
		-f -org -c eudetPPSFeb2014FEI4 -P:B_CheckRegion_MaskingRadius 1 | tee -a $LOGFILE
	
	echo "" | tee -a $LOGFILE

done

echo "script finished"
