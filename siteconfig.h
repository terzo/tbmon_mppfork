#include "tbconfig.h"
#include "stdlib.h"

//Site specific configuration. YOU WILL NEED TO EDIT THIS!
void trySet(char* config, char* path){
  if (strcmp(config, "notset") == 0 ){
    sprintf(config, "%s", path);
  }
  // change by Botho
  else {
    std::cout << "tryset: configuration was set already to " << config << std::endl << "value not changed" << endl;}
}

//Default output and extensions
#define SITECONFIG_SET
void siteConfig(TbConfig &config){
  //Path to where plots/logs will be written
//  trySet(config.outPath, (char*) "/remote/pcatlas62/testbeam/output/trunktest/");
//  trySet(config.outPath, (char*) "/remote/pcatlas62/testbeam/output/trunktest_feb/");
//  trySet(config.outPath, (char*) "/remote/pcatlas62/testbeam/output/new_reco_feb/trunk_nm2/");
//  trySet(config.outPath, (char*) "/remote/pcatlas70/testbeam/output/new_reco_feb/trunk20141020/");
  trySet(config.outPath, (char*) "/remote/pcatlas33/testbeam/new_analysis/");
//  trySet(config.outPath, (char*) "/mnt/scratch/testbeam/output/new_reco_feb/trunk20141021/");
//  trySet(config.outPath, (char*) "/remote/pcatlas70/testbeam/output/good_analysis/");
//Plot extension: .pdf, .png, .ps, .svg, or anything else recognized by ROOT
//set extension to none if no graphics files shall be generated
//  trySet(config.plotExtension, (char*) ".eps"); // commented out by Botho
    trySet(config.plotExtension, (char*) ".none");
//Default configuration
    trySet(config.tbslot, (char*) "eudetmar2013");
//ATLAS style plots (otherwise, use ROOT Plain style)
  // config.useAtlasStyle = true;
  // changed by Botho
  config.useAtlasStyle = false;
}

//default datapath for May 2008 BAT TB
#define SITEBAT2008_SET
void siteBAT2008(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2008/tbtrack");
}
//default datapath for May 2009 BAT TB
#define SITEBATMAY2009_SET
void siteBATMay2009(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2009-05/tbtrack.3");
}
#define SITEBATMAY2009OLD1_SET
void siteBATMay2009old1(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2009-05/tbtrack.1");
}
//default datapath for October 2009 BAT TB
#define SITEBATOCT2009_SET
void siteBATOct2009(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2009-10/tbtrack");
}
//default datapath for October 2009 EUDET TB
#define SITEEUDETOCT2009_SET
void siteEudetOct2009(TbConfig &config) {
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2009-10-eudet/tbtrack");
}
//default datapath for November 2009 EUDET TB
#define SITEEUDETNOV2009_SET
void siteEudetNov2009(TbConfig &config) {
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2009-11-eudet/tbtrack");
}
//Default datapath for 2010 simulation
#define SITESIMBAT2010_SET
void siteSimBat2010(TbConfig &config) {
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/sim");
  trySet(config.simDataPath, (char*) "/path/to/tbAnalysis/data/sim");
  config.isSimulation = true;
}
//default datapath for June 2010 EUDET TB
#define SITEEUDETJUNE2010_SET
void siteEudetJune2010(TbConfig &config) {
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2010-06-eudet/tbtrack");
}
//default datapath for Oct 2010 BAT TB
#define SITEBATOCT2010_SET
void siteBATOct2010(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2010-10/tbtrack");
}
//default datapath for November 2010 EUDET TB
#define SITEEUDETNOV2010_SET
void siteEudetNov2010(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2010-11-eudet/tbtrack");
}
//default datapath for February 2011 EUDET TB
#define SITEEUDETFEB2011_SET
void siteEudetFeb2011(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/tb2011-02-eudet/tbtrack");
}
//default datapath for June 2011 IBL EUDET TB
#define SITEEUDETJUN2011_SET
void siteEudetJun2011(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/cern_2011_jun_ibl/tbtrack");
}
//default datapath for July 2011 PPS EUDET TB
#define SITEEUDETPPSJUL2011_SET
void siteEudetPPSJul2011(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/cern_2011_jul_apix/tbtrack");
}
//default datapath for Sep 2011 IBL EUDET TB
#define SITEEUDETIBLSEP2011_SET
void siteEudetIBLSep2011(TbConfig &config){
  trySet(config.dataPath, (char*) "/afs/ipp/m/stte/eudet_data/cern_2011_sep_ibl/histo/");
}
//default datapath for Sep 2011 PPS EUDET TB
#define SITEEUDETPPSSEP2011_SET
void siteEudetPPSSep2011(TbConfig &config){
  trySet(config.dataPath, (char*) "/path/to/tbAnalysis/data/cern_2011_sep_apix/tbtrack");
}

//default datapath for September 2012 PPS EUDET TB at CERN
#define SITEEUDETPPSSEP2012FEI3_SET "/remote/pcatlas70/testbeam/data/cern_2012_sep_pps"
void siteEudetPPSSep2012FEI3(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSSEP2012FEI3_SET);
}

//default datapath for September 2012 PPS EUDET TB at CERN
#define SITEEUDETPPSSEP2012FEI4_SET "/afs/ipp-garching.mpg.de/home/s/stte/eudet_data/cern_2012_sep_pps/histo"
void siteEudetPPSSep2012FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSSEP2012FEI4_SET);
}

//default datapath for May 2012 PPS EUDET TB
#define SITEEUDETPPSMAY2012_SET "/remote/pcatlas33/testbeam/cern_2012_may_pps/ste_newreco"
void siteEudetPPSMay2012(TbConfig &config){
  //trySet(config.dataPath, (char*) "/remote/pcatlas33/testbeam/cern_2012_may_pps/results");
  //trySet(config.dataPath, (char*) "/remote/pcatlas33/testbeam/cern_2012_may_pps/stefano_rico");
  trySet(config.dataPath, (char*) SITEEUDETPPSMAY2012_SET);
}

//default datapath for Mar 2013 PPS EUDET TB
#define SITEEUDETPPSMAR2013FEI4_SET "/remote/pcatlas62/testbeam/data/desy_2013_mar_pps/histo/"
void siteEudetPPSMar2013FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSMAR2013FEI4_SET);
}

//default datapath for Aug 2013 PPS EUDET TB
#define SITEEUDETPPSAUG2013FEI34_SET "/remote/pcatlas62/testbeam/data_first_analysis/desy_2013_aug_pps/histo"
void siteEudetPPSAug2013FEI34(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSAUG2013FEI34_SET);
}

//default datapath for Nov 2013 PPS EUDET TB with FEI3 DUT (and FEI4 reference)
/* 
// for first reconstruction
#define SITEEUDETPPSNOV2013FEI34_SET "/remote/pcatlas62/testbeam/data_first_analysis/desy_2013_nov_pps/histo" */
// for more recently reconstruced 
#define SITEEUDETPPSNOV2013FEI34_SET "/remote/pcatlas70/testbeam/data/desy_2013_nov_pps/histo"
//#define SITEEUDETPPSNOV2013FEI34_SET "/mnt/scratch/testbeam/data/desy_2013_nov_pps/histo"
void siteEudetPPSNov2013FEI34(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSNOV2013FEI34_SET);
}

//default datapath for Nov 2013 PPS EUDET TB with FEI4 DUT (and FEI4 reference)
#define SITEEUDETPPSNOV2013FEI4_SET "/remote/pcatlas62/testbeam/data/desy_2013_nov_pps/histo"
//#define SITEEUDETPPSNOV2013FEI4_SET "/remote/pcatlas70/testbeam/data/desy_2013_nov_pps/histo"
void siteEudetPPSNov2013FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSNOV2013FEI4_SET);
}

//default datapath for Feb 2014 PPS EUDET TB with FEI4 DUT (and FEI4 reference)
//#define SITEEUDETPPSFEB2014FEI4_SET "/remote/pcatlas62/testbeam/data/desy_2014_feb_pps/histo"

#define SITEEUDETPPSFEB2014FEI4_SET "/mnt/scratch/testbeam/data/desy_2014_feb_pps/histo"
void siteEudetPPSFeb2014FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSFEB2014FEI4_SET);
}

//default datapath for Oct 2014 ITk EUDET TB with one FEI3 and two FEI4 DUT
#define SITEEUDETITKOCT2014FEI344_SET "/remote/pcatlas62/testbeam/data/cern_2014_oct_itk/histo"
void siteEudetITkOct2014FEI344(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETITKOCT2014FEI344_SET);
}

//default datapath for Mar 2015 ITk EUDET TB with one FEI4_25 and one FEI4 drtmund DUT
#define SITEEUDETITKMAR2015FEI4_25_SET "/remote/pcatlas62/testbeam/data/desy_2015_mar_itk/histo/"
void siteEudetITKMar2015FEI4_25(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETITKMAR2015FEI4_25_SET);
}

//default datapath for Mar 2015 ITk EUDET TB with one FEI4 and one FEI4 drtmund DUT
#define SITEEUDETITKMAR2015FEI4_SET "/remote/pcatlas62/testbeam/data/desy_2015_mar_itk/histo/"
void siteEudetITKMar2015FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETITKMAR2015FEI4_SET);
}

//default datapath for Jun 2012 IBL EUDET TB CERN
#define SITEEUDETIBLJUN2012FEI4_SET "/remote/pcatlas62/testbeam/data/cern_2012_jun_ibl/histo/"
void siteEudetIBLJun2012FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETIBLJUN2012FEI4_SET);
}

//default datapath for Aug 2013 PPS EUDET TB
#define SITEEUDETPPSAUG2013FEI4_SET "/afs/ipp-garching.mpg.de/home/s/stte/eudet_data/desy_2013_aug_pps/histo"
void siteEudetPPSAug2013FEI4(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETPPSAUG2013FEI4_SET);
}
//default datapath for Allpix
#define SITEEUDETALLPIX_SET "/remote/pcatlas2/testbeam/Simulations/ilcsoft/v01-17-03/Eutelescope/v00-09-01/jobsub/output/histo/"
void siteEudetAllpix(TbConfig &config){
  trySet(config.dataPath, (char*) SITEEUDETALLPIX_SET);
}
