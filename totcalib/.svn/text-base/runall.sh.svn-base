#!/bin/sh
option=$1
#if [ -z "$2" ]; then 
#	nev=99999999
#fi
if [ "${option}" = "all" ]
then
    ./totcal totcalresults/totcal_planar/2009-05-15-totclow_cal_0.tot 7 160 mask-PLANAR.txt 
    ./totcal totcalresults/totcal_planarnew/2009-07-14-Clo_0.tot 8.117 160new mask-PLANAR.txt 
    ./totcal totcalresults/totcal_sta/2009-05-15-totclow_cal_0.tot 7.1 164 mask-STA3D-3e.txt 
    ./totcal totcalresults/totcal_fbknew/2009-05-25-totclo_0.tot 7.581 168 mask-FBK3E7.txt 
    ./totcal totcalresults/totcal_fbknew/2009-05-25-totclo_0.tot 7.691 172 mask-FBK3EM5.txt 
    
#    ./tbmon runlistmay2009_bfieldoff_angle15.txt
fi

