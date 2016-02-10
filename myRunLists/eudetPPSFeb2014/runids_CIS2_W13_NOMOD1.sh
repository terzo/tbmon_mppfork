#! /bin/bash
# script needs to be called from trunk directory
# Botho Paschen
# 2014-04-01

# file for logging output of this script
LOGFILE=myRunLists/eudetPPSFeb2014/runids_CIS2_W13_NOMOD1.log



echo
echo "#############################################################"
echo "starting script" | tee $LOGFILE
echo
echo "output will be logged in "$LOGFILE
echo

for GEOID in 211 212 213 214 215 216 217 218
do

	echo "###################################################################" \
		| tee -a $LOGFILE
	echo "process geoid"$GEOID | tee -a $LOGFILE
	echo "###################################################################" \
		| tee -a $LOGFILE
	
	./bdrive/bdrive_args2.sh -l myRunLists/eudetPPSFeb2014/geoid$GEOID \
		eudetPPSFeb2014FEI4 | tee -a $LOGFILE
	
	echo "" | tee -a $LOGFILE

done

echo "script finished"
