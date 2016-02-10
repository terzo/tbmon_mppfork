#!/usr/bin/env python
"""
The goal of this program is to find a
good tuning for a simDut module by
varying its parameters, and comparing
the the resulting plots with real data.

Kyrre N. Sjobak (k.n.sjobak@fys.uio.no)
"""

usage = \
"""
Usage:
simdut_paramtuner -o <outName> -m <model> -P:param1 <value/range> -P:param2 <value/range> -c \"comment\"
"""

print "Welcome to simdut_paramtuner!"

class paramSet:
"""
@brief This class stores one single set of parameters, in the same format as the map TbConfig::cmdLineExtras<string,string>
    
It can also produce nicely formatted HTML tables of the current setup for use by the tuning report.
"""
    
    #Don't use static initializers
    # for non-static mutable member fields
    modelName = None
    argMap    = None 
    
    def __init__(self, modelName, argMap):
        self.modelName = modelName
        self.argMap    = argMap
    
    def makeHTML(self):
        """
        Produces an alfabetically sorted HTML table
        of the current parameter settings.
        """
        
        ret = """
        <table border=\"1\">
        <tr>
            <th> Parameter </th>
            <th> Setting </th>
        </tr>
        """
        for key in self.argMap.keys().sort():
            ret += """
            <tr>
                <td> %s </td>
                <td> %s </td> 
            </tr>
            """ % (key, self.argMap[key]);
        ret += """
        </table>
        """
        return ret;
    
    def __str__(self):
        return "Model = " + self.modelName + ", " + \
               "Parameters = " + str(self.argMap);  

#Configuration: Models and parameters (/w defaults)
modelParam_possible = {"Full3D_HP" : {"R_bias" : "7", "R_readout" : "7"} \
                        }

print "Possible models and their parameters:"
for (model,params) in modelParam_possible.iteritems():
    ps = paramSet(model, params)
    print "\t" + str(ps)

#Get tbmon setup by parsing siteconfig.h
siteconfig = open("siteconfig.h", "r")
siteconfig_lines = siteconfig.readlines();
siteconfig.close()

outPath = ""
rootStyle = ""
rootPalette = 0

for line in siteconfig_lines:
    if "config.outPath" in line:
        outPath = line.split("\"")[1] #Get whatever is between quotes
    elif "gROOT->SetStyle" in line:
        rootStyle = line.split("\"")[1]
    elif "gStyle->SetPalette" in line:
        #Get what is between "(" and ")", convert to int
        rootPalette = int(line.split("(")[1].split(")")[0])

del siteconfig_lines, siteconfig; 

print "Result from parsing siteConfig:"
print "\t outPath     = " + outPath
print "\t rootStyle   = " + rootStyle
print "\t rootPalette = " + str(rootPalette)

#Parse command line arguments
outname = ""
comment = ""
runs = "650 666"

print usage

#Get reference data histograms
refHistos = {}

#Setup HTML output

#Setup parameter sets
params = []

#Loop over parameter sets:
for param in params:
    pass

#Cleanup
