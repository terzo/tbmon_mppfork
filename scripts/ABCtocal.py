# -*- coding: utf-8 -*-



'''
Description: Small python script to convert TOT Calibs derived with USBPix for FE-I3 chips to the tbmon format.
Created: August 2011
Author: Philipp Weigell (pweigell@mpp.mpg.de)
Usage: Open a TOT_Calib Scan in DBeditor, save the 2d histograms for A, B & C as .C root macros (right click on the plot -> Save As)
       Copy files to the script dir, run sript, find the calibration in calib.out.
       So this for all duts you need, and cat the fiels together.
Caveat: Be sure to change line 22 of core/src/e4_totcalib.cc to handle UBPixlike fit constants instead of TurboDAQ fit constants an example would be:
       return ((tot/calibpp[index(col,row)]->calA)*calibpp[index(col,row)]->calC)-(calibpp[index(col,row)]->calB)/(1-(tot/calibpp[index(col,row)]->calA));

'''

dutnumber = 1 #Corresponds to tbmon dutnumber -10
a = open("A.C")
b = open("B.C")
c = open("C.C")

v_a ={}
v_b ={}
v_c ={}

for astr in a.readlines():
  if astr.find('SetBinContent')!=-1:
    pair = astr.split('(')[1].split(')')[0].split(',')
    v_a[(int(pair[0])%20,(int(pair[0])-int(pair[0])%20)/20)]=float(pair[1])
    
for bstr in b.readlines():
  if bstr.find('SetBinContent')!=-1:
    pair = bstr.split('(')[1].split(')')[0].split(',')
    #print bstr
    v_b[(int(pair[0])%20,(int(pair[0])-int(pair[0])%20)/20)]=float(pair[1])    

for cstr in c.readlines():
  if cstr.find('SetBinContent')!=-1:
    pair = cstr.split('(')[1].split(')')[0].split(',')
    #print cstr
    v_c[(int(pair[0])%20,(int(pair[0])-int(pair[0])%20)/20)]=float(pair[1])


out = open('calib.out','w')

for c in range(0,18):
  for r in range(0,160):
   try:
     print v_a[(c,r)],
   except KeyError:
     v_a[(c,r)]=0
   try:
     print v_b[(c,r)],
   except KeyError:
     v_b[(c,r)]=0
   try:
     print v_c[(c,r)]
   except KeyError:
     v_c[(c,r)]=0  
     
   out.write(str(dutnumber)+' '+str(c)+' '+str(r)+' 0 '+str(v_a[(c,r)])+' '+str(v_b[(c,r)])+' '+str(v_c[(c,r)])+"\n")
