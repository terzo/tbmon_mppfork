#include "driver.h"
#include "siteconfig.h"

//Event builders
#include "anglecuts.h"
#include "battrack.h"
#include "calcangles.h"
#include "checkregion.h"
#include "checkcentralregion.h"
#include "chi2builder.h"
#include "clusterdumper.h"
#include "clusterfinder.h"
#include "clustermasker.h"
#include "dutsync.h"
#include "batetacutter.h"
#include "eubuildtrack.h"
#include "eventbuilder.h"
#include "maskandlvl1.h"
#include "maskreader.h"
#include "lvl1cuts.h"
#include "pixelmasker.h"
#include "totcalibreader.h"
#include "translator.h"
#include "translatorRunningXY.h"
// sim stuff
#include "simBaseBuilder.h"
#include "simDutRunner.h"
#include "simPixelEdepBuilder.h"
#include "simTruthBuilder.h"

//Analysis headers
#include "batangledist.h"
#include "batunbiased.h"
#include "beamprofile.h"
#include "blank.h"
#include "botho.h"
#include "checkalign.h"
#include "checkalignRunningXY.h"
#include "batcheckdutsync.h"
#include "checktrack.h"
#include "clusterchecker.h"
#include "clustersvsrun.h"
#include "clusterpixtot.h"
#include "correlations.h"
#include "edgeefficiency.h"
#include "edgeefficiencyshift.h"
#include "efficiency.h"
#include "efficiency2.h"
#include "efficiencysimp.h"
#include "efficiencyvsrun.h"
#include "etawidth.h"
#include "getetacorr.h"
#include "hotpixelfinder.h"
#include "lvl1cut.h"
#include "maxcellresiduals.h"
#include "qEfficiency.h"
#include "qshare1D.h"
#include "qshare2D.h"
#include "readout.h"
#include "residuals.h"
#include "sumtot.h"
// sim stuff
#include "simDutEdep.h"
#include "simResiduals.h"
// graphics stuff
#include <TROOT.h>
#include <TStyle.h>
#include "AtlasStyle.h"

/*
 Check list for obtaining calibrations:
 -Run checkchi2 distribution, push chi2 builder on list with a reasonable cut
 -Run angledist over all runs with on DUT enabeled, push AngleCuts on list. 
 -Run hotpixelfinder with all DUTs enabeled over all runs, add masks to DUT
 -Run Check DUT sync with all DUTs over all runs, read back.
 -Run getetacorr with all DUTs enabeled oncs per angle/b-field config, read all eta calib files back in.
 -Run checktranslation over all runs for all good DUTs, pop bad runs from run list(investigate first)

 Notes:
 -Translatorevent builder should not be needed for newly reconstructed tracks.
 */

//Add all TbAnalyses objects, give it a proper and unique name so it is callable from the command line
void allAnalyses(TbConfig & config, DUT* dut) {
	//these guys produce text calib files
	//config.addAnalysis( new AngleDist, "angledist", dut );
	//config.addAnalysis( new CheckDUTSync, "checkdutsync", dut );
	//config.addAnalysis( new EtaWidth, "etawidth", dut);
	config.addAnalysis(new GetEtaCorr, "getetacorr", dut);
	config.addAnalysis(new HotPixelFinder, "hotpixelfinder", dut);
	//Regular analysis classes
	//config.addAnalysis( new batunbiased, "batunbiased", dut);
	config.addAnalysis(new BeamProfile, "beamprofile", dut);
	config.addAnalysis(new CheckAlign, "checkalign", dut);
	config.addAnalysis( new CheckTrack, "checktrack", dut);
	config.addAnalysis( new ClusterChecker, "clusterchecker", dut );
	config.addAnalysis( new ClustersVsRun, "clustervsrun", dut );
	config.addAnalysis(new ClusterPixToT, "clusterpixtot", dut);
	//config.addAnalysis( new Correlation, "correlations", dut);
	config.addAnalysis( new EdgeEfficiency, "edgeefficiency", dut );
	//config.addAnalysis( new EdgeEfficiencyShift, "edgeefficiencyshift", dut );
	config.addAnalysis(new Efficiency, "efficiency", dut);
	config.addAnalysis( new Efficiency2, "efficiency2", dut );
	//config.addAnalysis( new EfficiencySimp, "efficiencysimp", dut );
	config.addAnalysis(new EfficiencyVsRun, "efficiencyvsrun", dut);
	//config.addAnalysis( new Lvl1Cut, "lvl1cut", dut );
	config.addAnalysis(new MaxCellResiduals, "maxcellres", dut);
	config.addAnalysis(new qEfficiency, "qEfficiency", dut);
	config.addAnalysis(new QShare1D, "qshare1d", dut);
	config.addAnalysis(new QShare2D, "qshare2d", dut);
	//config.addAnalysis( new Readout, "readout" , dut );
	config.addAnalysis( new Residuals, "residuals" , dut );
	config.addAnalysis(new Botho, "botho", dut);
	config.addAnalysis(new SumTot, "sumtot", dut);
	//Simulation analysis requiring special simulation data
	//config.addAnalysis( new simResiduals, "simResiduals", dut);
	//config.addAnalysis( new simDutEdep, "simDutEdep", dut);
}

void EudetIBLSep2011(TbConfig &config) {
#ifdef SITEEUDETIBLSEP2011_SET
	siteEudetIBLSep2011(config);
#else
	cout << "Site config for EudetIBLSep2011 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	//config.MaxNumOfTrigToBeProcPerRun = 50;
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);
	eut->addMatchDUT(22);
	eut->addMatchDUT(23);

	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(30.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;	

	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	dutnumber = 22;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 23;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	

	dutnumber = 20;
	DUT* m20 = makeFEI4("VTT-NP2-20-E4", dutnumber, 2, 500.0, 25.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	
	dutnumber = 22;
	DUT* m22 = makeFEI4("VTT-NP2-22-E4", dutnumber, 2, 500.0, 25.0); 
	//DUT* m22 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m22->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m22->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m22->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m22->addToTCalib(calibfile.c_str());
	}
	config.addDut(m22);

	dutnumber = 23;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m23 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m23->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m23->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m23->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m23->addToTCalib(calibfile.c_str());
	}
	config.addDut(m23);

}

void EudetAllpix(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETALLPIX_SET
	siteEudetAllpix(config);
#else 
	cout << "Site config for EudetAllpix not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(0);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(20.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = false;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;	

	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	

	dutnumber = 20;
	DUT* m20 = makeFEI4("VTT-NP2-20-E4", dutnumber, 2, 500.0, 25.0, 10.0, 10.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}

void EudetIBLJun2012FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETIBLJUN2012FEI4_SET
	siteEudetIBLJun2012FEI4(config);
#else 
	cout << "Site config for EudetIBLJun2012FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);
	eut->addMatchDUT(22);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(20.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;	

	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	dutnumber = 22;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	

	dutnumber = 20;
	DUT* m20 = makeFEI4("VTT-NP2-20-E4", dutnumber, 2, 500.0, 25.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);
	
	dutnumber = 22;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m22 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m22->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m22->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m22->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m22->addToTCalib(calibfile.c_str());
	}
	config.addDut(m22);

/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}


void EudetITKMar2015FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETITKMAR2015FEI4_SET
	siteEudetITKMar2015FEI4(config);
#else 
	cout << "Site config for EudetITKMar2015FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(20.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;	

	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	

	dutnumber = 20;
	DUT* m20 = makeFEI4("VTT-NP2-20-E4", dutnumber, 2, 500.0, 25.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}

void EudetITKMar2015FEI4_25(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETITKMAR2015FEI4_25_SET
	siteEudetITKMar2015FEI4_25(config);
#else 
	cout << "Site config for EudetITKMar2015FEI4_25 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(20.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;	

	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, checkaligndir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	

	dutnumber = 20;
	DUT* m20 = makeFEI4_25("VTT-NP2-20-E4", dutnumber, 2, 500.0, 25.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}

void EudetPPSAug2013FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSAUG2013FEI4_SET
	siteEudetPPSAug2013FEI4(config);
#else 
	cout << "Site config for EudetPPSAug2013FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);
	eut->addMatchDUT(22);

	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0); 
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 600.0, 50.0);
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


void EudetPPSMar2013FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSMAR2013FEI4_SET
	siteEudetPPSMar2013FEI4(config);
#else 
	cout << "Site config for EudetPPSMar2013FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0); 
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 600.0, 50.0);
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


// October 2014 test beam at CERN SPS
// configuration with active edge FE-I3 and two FE-I4s
void EudetITkOct2014FEI344(TbConfig &config) {
#ifdef SITEEUDETITKOCT2014FEI344_SET
	siteEudetITkOct2014FEI344(config);
#else 
	cout << "Site config for EudetITkOct2014FEI344 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(20);
	eut->addMatchDUT(25);

	char outdir[800];
	strcpy(outdir, (const char*) config.outPath);
	

	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;
	
	dutnumber = 10;
	DUT* m10 = makeFEI3("FEI3", dutnumber, 2, 600.0, 50.0);
	m10->lv1Range(0, 16);
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 20;
	DUT* m20 = makeFEI4("VTT-NP2-20-E4", dutnumber, 2, 250.0, 50.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 25;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m25 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m25->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m25->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m25->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m25->addToTCalib(calibfile.c_str());
	}
	config.addDut(m25);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();
		
		dutnumber = 10;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 25;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void EudetPPSSep2012FEI3(TbConfig &config) {
#ifdef SITEEUDETPPSSEP2012FEI3_SET
	siteEudetPPSSep2012FEI3(config);
#else 
	cout << "Site config for EudetPPSSep2012FEI3 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);
	eut->addMatchDUT(12);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;
	DUT *tmpDUT;

	dutnumber = 10;
	DUT* m10 = makeFEI3("VTT-E3", dutnumber, 2, 600.0, 50.0);
	dutnumber = 11;
	DUT* m11 = makeFEI3("VTT-E5", dutnumber, 2, 600.0, 50.0);
	dutnumber = 12;
	DUT* m12 = makeFEI3("VTT-E6", dutnumber, 2, 600.0, 50.0);
	
	DUT *duts[3] = {m10, m11, m12};
	
	for(int i=0; i<3; i++) {
		tmpDUT = duts[i];
		dutnumber = 10+i;
		//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
		tmpDUT->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
		if (applyetacorr) {
			//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
			sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
					config.outPath, config.name, dutnumber);
			cout << "Applying eta corrections " << etacorrfile << " for DUT "
					<< dutnumber << endl;
			tmpDUT->ecorrs.addEcorr(etacorrfile);
		}
		if (masknoisyanddeadpixels) {
			//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
			sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
					config.name, dutnumber);
			cout << "Applying hot pixel mask " << maskfile << " for DUT "
					<< dutnumber << endl;
			tmpDUT->addMasks(maskfile);
		}
		if (addtotcalib) {
			calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
						<< dutnumber << endl;
			tmpDUT->addToTCalib(calibfile.c_str());
		}
		config.addDut(tmpDUT);
	}

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		for(int i=0; i<3; i++) {
			dutnumber = 10+i;
			sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
				config.outPath, checkaligndir, config.name, dutnumber);
			cout << "Applying translations " << alignfile << " for DUT "
					<< dutnumber << endl;
			translator->addTranslation(alignfile, dutnumber);
		}

		config.addBuilder(translator);
	}
}

void EudetPPSSep2012FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSSEP2012FEI4_SET
	siteEudetPPSSep2012FEI4(config);
#else 
	cout << "Site config for EudetPPSSep2012FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = false;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = false;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0); 
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 600.0, 50.0);
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


void EudetPPSFeb2014FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSFEB2014FEI4_SET
	siteEudetPPSFeb2014FEI4(config);
#else 
	cout << "Site config for EudetPPSFeb2014FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	DUT* m20 = makeFEI4("VTT-NP2-20-E4", dutnumber, 2, 250.0, 50.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("CIS2-W02-IBL3", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSNov2013FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSNOV2013FEI4_SET
	siteEudetPPSNov2013FEI4(config);
#else 
	cout << "Site config for EudetPPSNov2013FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	DUT* m20 = makeFEI4("CI2-W02-IBL2", dutnumber, 2, 250.0, 50.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m21 = makeFEI4("VTT-NP1-8-E4", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSNov2013FEI34(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSNOV2013FEI34_SET
	siteEudetPPSNov2013FEI34(config);
#else 
	cout << "Site config for EudetPPSNov2013FEI34 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(20);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("VTT-NP1-4-E5", dutnumber, 2, 600.0, 50.0); 
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 336.0, 80.0);
	m10->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 20;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 400.0, 50.0);
	DUT* m20 = makeFEI4("CIS-67-W2-IBL2", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void CBJuly2013(TbConfig &config) {
#ifdef SITECBJULY2013_SET
	siteCBJuly2013(config);
#else 
	cout << "Site config for CBJuly2013 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);

	eut->nMatches(0);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibdir = config.outPath;
	std::string calibfile = calibdir + "cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	DUT* m20 = makeFEI4("SCC-167", dutnumber, 2, 250.0, 50.0); // makeFEI4(DUT_Name, DUT_ID, 2, edge_pixel_size_Xum, edge_pixel_size_Yum)
	m20->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		if ((config.runList.front() >= 61249)
				&& (config.runList.back() <= 61275)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61274_geoID40.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61276)
				&& (config.runList.back() <= 61289)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID42.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61290)
				&& (config.runList.back() <= 61297)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID43.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61298)
				&& (config.runList.back() <= 61303)) {
			//should have an own calibration file - missing
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61304)
				&& (config.runList.back() <= 61314)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		/*
		 // dut13, average shift from checkalign analysis
		 translator->addTranslation(13, 30003, 30087, -28.639, -13.53);
		 translator->addTranslation(13, 30724, 30755, -13.276, -8.044);
		 translator->addTranslation(13, 30759, 30813, -12.7991, 2.1955);
		 translator->addTranslation(13, 30815, 30856, -12.2124, 3.2513);
		 translator->addTranslation(13, 30864, 30909, -13.0827, 4.3955);
		 translator->addTranslation(13, 30936, 30989, -13.2130, 3.2718);
		 translator->addTranslation(13, 30995, 31032, -28.6005, 1.7431);
		 // dut14, average shift from checkalign analysis
		 translator->addTranslation(14, 30003, 30087, -16.2673, 5.0409);
		 translator->addTranslation(14, 30724, 30755, -4.2859, -0.1675);
		 translator->addTranslation(14, 30759, 30813, -3.9610, -1.6993);
		 translator->addTranslation(14, 30815, 30856, -1.7750, -7.5424);
		 translator->addTranslation(14, 30864, 30909, -2.6596, -8.2410);
		 translator->addTranslation(14, 30936, 30989, -7.5564, -2.0348);
		 translator->addTranslation(14, 30995, 31032, -14.4534, -0.4106);
		 */

		config.addBuilder(translator);
	}
}

void EudetPPSSep2011(TbConfig &config) {
#ifdef SITEEUDETPPSSEP2011_SET
	siteEudetPPSSep2011(config);
#else
	cout << "Site config for EudetPPSSep2011 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);
	eut->addMatchDUT(12);
	eut->addMatchDUT(13);

	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = false;
	// apply eta corrections?
	bool applyetacorr = false;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = false;

	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibdir = "PPScalibs";
	std::string calibfile = calibdir + "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 400.0, 50.0); // makeFEI4(DUT_Name, DUT_ID, 2, edge_pixel_size_Xum, edge_pixel_size_Yum)
	m10->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		if ((config.runList.front() >= 61249)
				&& (config.runList.back() <= 61275)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61274_geoID40.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61276)
				&& (config.runList.back() <= 61289)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID42.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61290)
				&& (config.runList.back() <= 61297)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID43.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61298)
				&& (config.runList.back() <= 61303)) {
			//should have an own calibration file - missing
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61304)
				&& (config.runList.back() <= 61314)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 11;
	DUT* m11 = makeFEI3("DO-I-11", dutnumber, 2, 400.0, 50.0);
	m11->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m11->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m11->addMasks(maskfile);
	}
	if (addtotcalib) {
		if ((config.runList.front() >= 61249)
				&& (config.runList.back() <= 61275)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61274_geoID40.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61276)
				&& (config.runList.back() <= 61289)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID42.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61290)
				&& (config.runList.back() <= 61297)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID43.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61298)
				&& (config.runList.back() <= 61303)) {
			//should have an own calibration file - missing
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61304)
				&& (config.runList.back() <= 61314)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		m11->addToTCalib(calibfile.c_str());
	}
	config.addDut(m11);

	dutnumber = 12;
	DUT* m12 = makeFEI3("DO-I-5", dutnumber, 2, 400.0, 50.0);
	m12->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m12->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m12->addMasks(maskfile);
	}
	if (addtotcalib) {
		if ((config.runList.front() >= 61249)
				&& (config.runList.back() <= 61275)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61274_geoID40.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61276)
				&& (config.runList.back() <= 61289)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID42.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61290)
				&& (config.runList.back() <= 61297)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID43.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61298)
				&& (config.runList.back() <= 61303)) {
			//should have an own calibration file - missing
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61304)
				&& (config.runList.back() <= 61314)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		m12->addToTCalib(calibfile.c_str());
	}
	config.addDut(m12);

	dutnumber = 13;
	DUT* m13 = makeFEI3("DO-I-12", dutnumber, 2, 400.0, 50.0);
	m13->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m13->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m13->addMasks(maskfile);
	}
	if (addtotcalib) {
		if ((config.runList.front() >= 61249)
				&& (config.runList.back() <= 61275)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61274_geoID40.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61276)
				&& (config.runList.back() <= 61289)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID42.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61290)
				&& (config.runList.back() <= 61297)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61289_geoID43.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61298)
				&& (config.runList.back() <= 61303)) {
			//should have an own calibration file - missing
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		if ((config.runList.front() >= 61304)
				&& (config.runList.back() <= 61314)) {
			calibfile = "PPScalibs/totcalibs/scan_after_run61314_geoID45.out";
			cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		};
		m13->addToTCalib(calibfile.c_str());
	}
	config.addDut(m13);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 11;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 12;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 13;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		/*
		 // dut13, average shift from checkalign analysis
		 translator->addTranslation(13, 30003, 30087, -28.639, -13.53);
		 translator->addTranslation(13, 30724, 30755, -13.276, -8.044);
		 translator->addTranslation(13, 30759, 30813, -12.7991, 2.1955);
		 translator->addTranslation(13, 30815, 30856, -12.2124, 3.2513);
		 translator->addTranslation(13, 30864, 30909, -13.0827, 4.3955);
		 translator->addTranslation(13, 30936, 30989, -13.2130, 3.2718);
		 translator->addTranslation(13, 30995, 31032, -28.6005, 1.7431);
		 // dut14, average shift from checkalign analysis
		 translator->addTranslation(14, 30003, 30087, -16.2673, 5.0409);
		 translator->addTranslation(14, 30724, 30755, -4.2859, -0.1675);
		 translator->addTranslation(14, 30759, 30813, -3.9610, -1.6993);
		 translator->addTranslation(14, 30815, 30856, -1.7750, -7.5424);
		 translator->addTranslation(14, 30864, 30909, -2.6596, -8.2410);
		 translator->addTranslation(14, 30936, 30989, -7.5564, -2.0348);
		 translator->addTranslation(14, 30995, 31032, -14.4534, -0.4106);
		 */

		config.addBuilder(translator);
	}
}


void EudetPPSJuly2011(TbConfig &config) {
#ifdef SITEEUDETPPSJULY2011_SET
	siteEudetPPSJuly2011(config);
#else 
	cout << "Site config for EudetPPSJuly2011 not found" << endl;
#endif
	cout << "Configuration not yet available" << endl;
}

void EudetJune2011(TbConfig &config) {
#ifdef SITEEUDETJUNE2011_SET
	siteEudetJune2011(config);
#else
	cout << "Site config for EudetJune2011 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);
	eut->addMatchDUT(22);
	eut->addMatchDUT(23);

	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new LVL1Cuts);
	//config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new ClusterMasker);
	//  config.addBuilder(new CheckRegion);
	config.addBuilder(new CheckCentralRegion);

	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	// bool masknoisyanddeadpixels = false;
	bool masknoisyanddeadpixels = true;

	// apply eta corrections?
	bool applyetacorr = false;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = false;

	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibdir = "IBLcalibs";
	std::string calibfile = calibdir + "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	DUT* m20 = makeFEI4("SCC-41", dutnumber, 2, 250.0, 50.0); // makeFEI4(DUT_Name, DUT_ID, 2, edge_pixel_size_Xum, edge_pixel_size_Yum)
	m20->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);

		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying masked pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		//  m20->addMasks(maskfile);
		m20->addMaskedPixels(maskfile, 1);

	}
	if (addtotcalib) {
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
				<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCC-60", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying masked pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		//    m21->addMasks(maskfile);
		m21->addMaskedPixels(maskfile, 1);

	}
	if (addtotcalib) {
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
				<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	dutnumber = 22;
	DUT* m22 = makeFEI4("SCC-34", dutnumber, 2, 250.0, 50.0);
	m22->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m22->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying masked pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		//    m22->addMasks(maskfile);
		m22->addMaskedPixels(maskfile, 1);

	}
	if (addtotcalib) {
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
				<< dutnumber << endl;
		m22->addToTCalib(calibfile.c_str());
	}
	config.addDut(m22);

	dutnumber = 23;
	DUT* m23 = makeFEI4("SCC-40", dutnumber, 2, 250.0, 50.0);
	m23->lv1Range(0, 10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%s/%s-getetacorr-%d-etacalib.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m23->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%s/%s-hotpixelfinder-%d-masks.txt", calibdir.c_str(),
				config.name, dutnumber);
		cout << "Applying masked pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		//    m23->addMasks(maskfile);
		m23->addMaskedPixels(maskfile, 1);

	}
	if (addtotcalib) {
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
				<< dutnumber << endl;
		m23->addToTCalib(calibfile.c_str());
	}
	config.addDut(m23);

	/*
	 dutnumber = 24;
	 DUT* m24 = makeFEI4("SCC-55", dutnumber, 2);
	 m24->lv1Range(0,10); // these cuts are only needed to look for hot pixels (in out-of-time events)
	 if (applyetacorr) {
	 // adding eta correction file. You may get it running the analysis getetacorr
	 //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
	 sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",calibdir.c_str(),config.name,dutnumber);
	 cout << "Applying eta corrections " << etacorrfile << " for DUT " << dutnumber << endl;
	 m24->ecorrs.addEcorr(etacorrfile);
	 }
	 if (masknoisyanddeadpixels) {
	 // adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
	 //sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
	 sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",calibdir.c_str(),config.name,dutnumber);
	 cout << "Applying hot pixel mask " << maskfile << " for DUT " << dutnumber << endl;
	 m24->addMasks(maskfile);
	 }
	 if (addtotcalib) {
	 cout << "Applying TOT calib from file " << calibfile << " for DUT " << dutnumber << endl;
	 m24->addToTCalib(calibfile.c_str());
	 }
	 config.addDut(m24);
	 */

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 13;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 14;
		sprintf(alignfile, "%s/%s-checkalign-%d-translations.txt",
				calibdir.c_str(), config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		/*
		 // dut13, average shift from checkalign analysis
		 translator->addTranslation(13, 30003, 30087, -28.639, -13.53);
		 translator->addTranslation(13, 30724, 30755, -13.276, -8.044);
		 translator->addTranslation(13, 30759, 30813, -12.7991, 2.1955);
		 translator->addTranslation(13, 30815, 30856, -12.2124, 3.2513);
		 translator->addTranslation(13, 30864, 30909, -13.0827, 4.3955);
		 translator->addTranslation(13, 30936, 30989, -13.2130, 3.2718);
		 translator->addTranslation(13, 30995, 31032, -28.6005, 1.7431);
		 // dut14, average shift from checkalign analysis
		 translator->addTranslation(14, 30003, 30087, -16.2673, 5.0409);
		 translator->addTranslation(14, 30724, 30755, -4.2859, -0.1675);
		 translator->addTranslation(14, 30759, 30813, -3.9610, -1.6993);
		 translator->addTranslation(14, 30815, 30856, -1.7750, -7.5424);
		 translator->addTranslation(14, 30864, 30909, -2.6596, -8.2410);
		 translator->addTranslation(14, 30936, 30989, -7.5564, -2.0348);
		 translator->addTranslation(14, 30995, 31032, -14.4534, -0.4106);
		 */

		config.addBuilder(translator);
	}
}

void EudetFeb2011(TbConfig &config) {
#ifdef SITEEUDETFEB2011_SET
	siteEudetFeb2011(config);
#else
	cout << "Site config for EudetFeb2011 not found" << endl;
#endif
	cout << "Configuration not yet available" << endl;
}

void EudetNov2010(TbConfig &config) {
#ifdef SITEEUDETNOV2010_SET
	siteEudetNov2010(config);
#else 
	cout << "Site config for EudetNov2010 not found" << endl;
#endif
	cout << "dataPath: " << config.dataPath << endl;
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);
	eut->addMatchDUT(12);
	//eut->addMatchDUT(13);
	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	//DUT* m10 = new DUT("STA3D", 10, 3); // BJD
	DUT* m10 = makeFEI3("STA3D", 10, 3); // makeFEI3(DUT_Name, DUT_ID, 3, edge_pixel_size_Xum, edge_pixel_size_Yum)
	m10->lv1Range(1, 10);
	//m10->addMasks("calibs/*.txt");
	//m10->ecorrs.addEcorr("calibs/*.txt");
	config.addDut(m10);
	//DUT* m11 = new DUT("PPS1", 11, 3); // BJD
	DUT* m11 = makeFEI3("PPS1", 11, 3);
	m11->lv1Range(1, 10);
	//m11->addMasks("calibs/*.txt");
	//m11->ecorrs.addEcorr("calibs/*.txt");
	config.addDut(m11);
	//DUT* m12 = new DUT("PPS0", 12, 3);
	DUT* m12 = makeFEI3("PPS0", 12, 3);
	m12->lv1Range(1, 10);
	//m12->addMasks("calibs/*.txt");
	//m12->ecorrs.addEcorr("calibs/*.txt");
	config.addDut(m12);
}

void EudetJuly2010(TbConfig& config) {
#ifdef SITEEUDETJULY2010_SET
	siteEUDETJuly2010(config);
#else 
	cout << "Site config for EudetJuly2010 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();

	//
	// add here all devices in your run
	//
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);
	eut->addMatchDUT(12);
	eut->addMatchDUT(13);
	eut->addMatchDUT(14);
	eut->addMatchDUT(15);
	eut->addMatchDUT(16);
	eut->addMatchDUT(17);

	//eut->addTranslation(13, 7200,7210,-340,-420);
	eut->nMatches(1);

	config.addBuilder(eut);
	//config.addBuilder(new ToTCalibReader);
	//config.addBuilder(new MaskReader);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder());
	config.addBuilder(new CheckRegion);
	//config.addBuilder(new MaskAndLvl1);
	//config.addBuilder(new ClusterDumper());

	//Translator* translator = new Translator();

	//eut->addTranslation(12, 0,100000,-130,00);
	//config.eventbuilders.push_back(translator);

	//
	// initialize each dut
	//
	int run = 11193;
	std::string calibdir = "calibs/";
	std::string calibfile =
			calibdir
					+ "ToT_CLOW_hv_2010-07-14-BCN-DO6-DO3-MPP1-D08-MPP2-DO9-DO7-3200e-60tot-20ke_DD89_800V_cal.out";
	cout << "ToT calibration file = " << calibfile.c_str() << endl;

	bool addtotcalib = true;
	bool applyetacorr = true;
	char etacorrfile[200];
	bool masknoisyanddeadpixels = true;
	char maskfile[200];

	DUT* m10 = makeFEI3("BCN", 10, 3); // constructor: DUT(name, ID, 3)
	m10->lv1Range(0, 10); // cuts for LV1 distribution

	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 10);
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				10);
		m10->addMasks(maskfile);
	}
	if (addtotcalib)
		m10->addToTCalib(calibfile.c_str());
	config.addDut(m10);

	DUT* m11 = makeFEI3("DO6", 11, 3);
	m11->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 11);
		m11->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				11);
		m11->addMasks(maskfile);
	}
	if (addtotcalib)
		m11->addToTCalib(calibfile.c_str());
	config.addDut(m11);

	DUT* m12 = makeFEI3("DO3", 12, 3);
	m12->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 12);
		m12->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				12);
		m12->addMasks(maskfile);
	}
	if (addtotcalib)
		m12->addToTCalib(calibfile.c_str());
	config.addDut(m12);

	DUT* m13 = makeFEI3("MPP1", 13, 3);
	m13->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 13);
		m13->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				13);
		m13->addMasks(maskfile);
	}
	if (addtotcalib)
		m13->addToTCalib(calibfile.c_str());
	config.addDut(m13);

	DUT* m14 = makeFEI3("DO8", 14, 3);
	m14->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 14);
		m14->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				14);
		m14->addMasks(maskfile);
	}
	if (addtotcalib)
		m14->addToTCalib(calibfile.c_str());
	config.addDut(m14);

	DUT* m15 = makeFEI3("MPP2", 15, 3);
	m15->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 15);
		m15->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				15);
		m15->addMasks(maskfile);
	}
	if (addtotcalib)
		m15->addToTCalib(calibfile.c_str());
	config.addDut(m15);

	DUT* m16 = makeFEI3("DO9", 16, 3);
	m16->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 16);
		m16->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				16);
		m16->addMasks(maskfile);
	}
	if (addtotcalib)
		m16->addToTCalib(calibfile.c_str());
	config.addDut(m16);

	DUT* m17 = makeFEI3("DO7", 17, 3);
	m17->lv1Range(0, 10);
	if (applyetacorr) {
		sprintf(etacorrfile, "CERN10_out/run%06d-getetacorr-%d-etacalib.txt",
				run, 17);
		m17->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		sprintf(maskfile, "CERN10_out/run%06d-hotpixelfinder-%d-masks.txt", run,
				17);
		m17->addMasks(maskfile);
	}
	if (addtotcalib)
		m17->addToTCalib(calibfile.c_str());
	config.addDut(m17);

	/*
	 DUT* m1 = new DUT("TEL", 1, 0);
	 config.addDut(m1);
	 DUT* m2 = new DUT("TEL", 2, 0);
	 config.addDut(m2);
	 DUT* m3 = new DUT("TEL", 3, 0);
	 config.addDut(m3);
	 DUT* m4 = new DUT("TEL", 4, 0);
	 config.addDut(m4);
	 DUT* m5 = new DUT("TEL", 5, 0);
	 config.addDut(m5);
	 */

	//Translator* translator = new Translator();
	//translator->addTranslation(16,0, 9999, 23.12,-73.58);
	//14
	//translator->addTranslation(14, 8220, 8220, -1500, 0);
	//15
	//translator->addTranslation(15, 8220, 8220,-2000,-200);
	//config.eventbuilders.push_back(translator);
}

void EudetJune2010(TbConfig &config) {
#ifdef SITEEUDETJUNE2010_SET
	siteEudetJune2010(config);
#else 
	cout << "Site config for EudetJune2010 not found" << endl;
#endif
	cout << config.dataPath << endl;
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);
	eut->addMatchDUT(12);
	//eut->addMatchDUT(13);
	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	DUT* m10 = new DUT("PLA", 10, 3);
	m10->lv1Range(1, 9);
	m10->addMasks("calibs/june2010-eudet-all-hotpixelfinder-10-masks.txt");
	m10->ecorrs.addEcorr(
			"calibs/june2010-eudet-0degree-getetacorr-10-etacalib.txt");
	m10->ecorrs.addEcorr(
			"calibs/june2010-eudet-15degree-getetacorr-10-etacalib.txt");
	config.addDut(m10);

	DUT* m11 = new DUT("FBK", 11, 3);
	m11->lv1Range(0, 10);
	m11->addMasks("calibs/june2010-eudet-all-hotpixelfinder-11-masks.txt");
	m11->ecorrs.addEcorr(
			"calibs/june2010-eudet-0degree-getetacorr-11-etacalib.txt");
	m11->ecorrs.addEcorr(
			"calibs/june2010-eudet-15degree-getetacorr-11-etacalib.txt");
	config.addDut(m11);

	DUT* m12 = new DUT("FBK", 12, 3);
	m12->lv1Range(0, 10);
	m12->addMasks("calibs/june2010-eudet-all-hotpixelfinder-12-masks.txt");
	m12->ecorrs.addEcorr(
			"calibs/june2010-eudet-0degree-getetacorr-12-etacalib.txt");
	m12->ecorrs.addEcorr(
			"calibs/june2010-eudet-15degree-getetacorr-12-etacalib.txt");
	config.addDut(m12);

	DUT* m13 = new DUT("FBK", 13, 3);
	m13->lv1Range(2, 14);
	//m13->addMasks("calibs/dut13-masks.txt");
	//m13->ecorrs.addEcorr("calibs/june2010_eudet_0degree-getetacorr-13-etacalib.txt");
	//m13->ecorrs.addEcorr("calibs/june2010_eudet_15degree-getetacorr-13-etacalib.txt");
	config.addDut(m13);
}

void EudetPPSMay2012FEI4(TbConfig &config) {
  //Used for geoID's 50-55 where both modules are SCM
#ifdef SITEEUDETPPSMAY2012FEI4_SET
	siteEudetPPSMay2012FEI4(config);
#else 
	cout << "Site config for EudetPPSMay2012FEI4 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(20);
	eut->addMatchDUT(21);


	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = true;
	// apply eta corrections?
	bool applyetacorr = true;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = true;

	//true if ./setup was run, false if not
	bool secondtranslation = true;
	int checkaligndir;
	if (secondtranslation) {
	  checkaligndir=2;
	}
	else {
	  checkaligndir=1;
	}
	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	//std::string currentdir=config.name;
	//std::string calibdir = "/mnt/scratch/mellenburg/"+currentdir;
	std::string calibfile = "/cal.out";
	//cout << "ToT calibration file = " << calibfile.c_str() << endl;

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0); 
	DUT* m20 = makeFEI4("SCM127", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(3, 9); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		//sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		//sprintf(maskfile,"%s/%s-hotpixelfinder-%d-masks.txt",config.outPath,config.name,dutnumber);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	//DUT* m11 = makeFEI3("SLID9", dutnumber, 2, 600.0, 50.0);
	DUT* m21 = makeFEI4("SCM130", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(3, 9);
	if (applyetacorr) {
	  //sprintf(etacorrfile,"%s/%s-getetacorr-%d-etacalib.txt",config.outPath,config.name,dutnumber);
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				config.outPath, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
	  //sprintf(maskfile,"%shotpixelfinder/hotpix+deadpix-mask-11.txt",config.outPath);
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", config.outPath, config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		
		calibfile="/mnt/scratch/mellenburg/othercalib/SLID9.out";
		cout << "Applying TOT calib from file " << calibfile << " for DUT "
					<< dutnumber << endl;
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
			config.outPath, checkaligndir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetNov2009(TbConfig& config) {
#ifdef SITEEUDETNOV2009_SET
	siteEudetNov2009(config);
#else 
	cout << "Site config for EudetNov2009 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(14);
	eut->addMatchDUT(15);
	eut->addMatchDUT(17);
	eut->nMatches(0);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	DUT* m14 = new DUT("PLA", 14, 3);
	m14->lv1Range(0, 10);
	config.addDut(m14);
	DUT* m15 = new DUT("STA", 15, 3);
	m15->lv1Range(0, 10);
	config.addDut(m15);
	DUT* m17 = new DUT("SIN", 17, 3);
	m17->lv1Range(0, 10);
	config.addDut(m17);
}

void EudetOct2009(TbConfig& config) {
#ifdef SITEEUDETOCT2009_SET
	siteEudetOct2009(config);
#else 
	cout << "Site config for EudetOct2009 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	config.addBuilder(new CheckRegion);
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);
	eut->addMatchDUT(13);
	eut->nMatches(1);
	void addTranslation(int iden, int firstRun, int lastRun, double shiftX,
			double shiftY);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	DUT* m10 = new DUT("FBK", 10, 3);
	config.addDut(m10);
	DUT* m11 = new DUT("STA", 11, 3);
	config.addDut(m11);
	DUT* m13 = new DUT("PLA", 13, 3);
	config.addDut(m13);
	DUT* m1 = new DUT("TEL", 1, 0);
	config.addDut(m1);
	DUT* m2 = new DUT("TEL", 2, 0);
	config.addDut(m2);
	DUT* m3 = new DUT("TEL", 3, 0);
	config.addDut(m3);
	DUT* m4 = new DUT("TEL", 4, 0);
	config.addDut(m4);
	DUT* m5 = new DUT("TEL", 5, 0);
	config.addDut(m5);
}

void BAT2008(TbConfig &config) {
#ifdef SITEBAT2008_SET
	siteBAT2008(config);
#else 
	cout << "Site config for Bat2008 not found" << endl;
#endif
	config.addBuilder(new BatTrack);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder());
	config.addBuilder(new CheckRegion);

	Chi2Builder* chi2 = new Chi2Builder(15.0); //Chi2 cut at 20.
	config.addBuilder((EventBuilder*) chi2);
	DUT* m160 = new DUT("PLANAR", 160, 0);
	//m162->lv1Range(3,10);
	//m162->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-162-masks.txt");
	//m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-162-etacalib.txt");
	config.addDut(m160);

	DUT* m164 = new DUT("STA-3G", 164, 0);
	//m163->lv1Range(3,10);
	//m163->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-163-masks.txt");
	//m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-163-etacalib.txt");
	config.addDut(m164);

}

void may2009DUTs(TbConfig &config) {
	//Configure DUTs used in the may2009 test beam
	DUT* m160 = new DUT("PLANAR", 160, 0);
	//Range of acceptable lv1 values, keep loose to avoid false positives in noise scan
	m160->lv1Range(2, 8);
	m160->addMasks("calibs/may2009-hotpixelfinder-160-masks.txt");
	//Eta corrections, output from getEtaCorr.
	m160->ecorrs.addEcorr(
			"calibs/may2009-boff-a15-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-1-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-2-getetacorr-160-etacalib.txt");
	m160->addToTCalib("calibs/totToQ_160.root");
	//FIXME
	//m160->totcalib->chargeCorrection(8.117 / 7.0); //wrong C-low
	config.addDut(m160);
	DUT* m164 = new DUT("STA-3G", 164, 0);
	m164->lv1Range(2, 10);
	m164->addMasks("calibs/may2009-hotpixelfinder-164-masks.txt");
	m164->ecorrs.addEcorr(
			"calibs/may2009-boff-a15-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-1-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-2-getetacorr-164-etacalib.txt");
	m164->addToTCalib("calibs/totToQ_164.root");
	config.addDut(m164);
	DUT* m168 = new DUT("FBK-3E7", 168, 0);
	m168->lv1Range(3, 11);
	m168->addMasks("calibs/may2009-hotpixelfinder-168-masks.txt");
	m168->ecorrs.addEcorr(
			"calibs/may2009-boff-a15-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-1-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-2-getetacorr-168-etacalib.txt");
	m168->addToTCalib("calibs/totToQ_168.root");
	config.addDut(m168);
	DUT* m172 = new DUT("FBK-3EM5", 172, 0);
	m172->lv1Range(2, 10);
	m172->addMasks("calibs/may2009-hotpixelfinder-172-masks.txt");
	m172->ecorrs.addEcorr(
			"calibs/may2009-boff-a15-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-1-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr(
			"calibs/may2009-boff-a0-2-getetacorr-172-etacalib.txt");
	m172->addToTCalib("calibs/totToQ_172.root");
	config.addDut(m172);
}

void may2009DUTsOld1(TbConfig &config) {
	//Configure DUTs used in the may2009 test beam
	DUT* m160 = new DUT("PLANAR", 160, 0);
	m160->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-160-etacalib.txt");
	m160->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-160-etacalib.txt");
	m160->addToTCalib("calibs/totToQ_160.root");
	//FIXME:
	//m160->totcalib->chargeCorrection(8.117 / 7.0); //wrong C-low
	config.addDut(m160);
	DUT* m164 = new DUT("STA-3G", 164, 0);
	m164->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-164-etacalib.txt");
	m164->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-164-etacalib.txt");
	m164->addToTCalib("calibs/totToQ_164.root");
	config.addDut(m164);
	DUT* m168 = new DUT("FBK-3E7", 168, 0);
	m168->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-168-etacalib.txt");
	m168->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-168-etacalib.txt");
	m168->addToTCalib("calibs/totToQ_168.root");
	config.addDut(m168);
	DUT* m172 = new DUT("FBK-3EM5", 172, 0);
	m172->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-172-etacalib.txt");
	m172->ecorrs.addEcorr(
			"calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-172-etacalib.txt");
	m172->addToTCalib("calibs/totToQ_172.root");
	config.addDut(m172);
}

//Configuration for BAT May2009 
void BATmay2009(TbConfig &config) {
	cout << "Using configuration may2009" << endl;
	//Site configuration
#ifdef SITEBATMAY2009_SET
	siteBATMay2009(config);
#else
	cout << "Site config for BATmay2009 not found" << endl;
#endif
	// Set up and configure eventbuilders
	config.addBuilder(new BatTrack);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	Chi2Builder* chi2 = new Chi2Builder(15.0); //Chi2 cut
	config.addBuilder((EventBuilder*) chi2);
	AngleCuts* angles = new AngleCuts("calibs/may2009-angledist.txt", 1.5); //Angle cut at mean +/- 1.5 sigma
	config.eventbuilders.push_back((EventBuilder*) angles);
	DutSync* dutsync = new DutSync();
	dutsync->readFile("calibs/dutsync-checkdutsync-160-dutsync.txt");
	config.addBuilder(dutsync);
	may2009DUTs(config);
}

//Configuration for BAT May2009 tbtrack.1 files
void BATmay2009old1(TbConfig &config) {
	cout << "Using configuration may2009.tbtrack.1" << endl;
	//Site configuration
#ifdef SITEBATMAY2009OLD1_SET
	siteBATMay2009old1(config);
#else
	cout << "Site config for BATmay2009old1 not found" << endl;
#endif
	// Set up and configure eventbuilders
	config.addBuilder(new BatTrack);
	config.addBuilder(new PixelMasker); // need this to avoid empty hits vector!
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	Chi2Builder* chi2 = new Chi2Builder(20.0); //Chi2 cut
	config.addBuilder((EventBuilder*) chi2);
	AngleCuts* angles = new AngleCuts(
			"calibs/unprocessed/may2009-unprocessed-angledist.txt", 1.5); //Angle cut at mean +/- 1.5 sigma
	config.addBuilder((EventBuilder*) angles);
	DutSync* dutsync = new DutSync();
	dutsync->readFile(
			"calibs/unprocessed/run500-1500-checkdutsync-160-dutsync.txt");
	config.addBuilder(dutsync);
	Translator* translator = new Translator();
	//160
	translator->addTranslation(160, 600, 803, 6.0, 1.23e-1);
	translator->addTranslation(160, 804, 998, 12.7, -1.69);
	translator->addTranslation(160, 999, 1204, 6.85, 5.03e-1);
	translator->addTranslation(160, 1205, 1307, -3.01, -5.55e-2);
	translator->addTranslation(160, 1308, 1364, -1.79, 2.11e-1);
	//164
	translator->addTranslation(164, 600, 803, -6.9, 3.36e-1);
	translator->addTranslation(164, 804, 998, 13.1, -1.32e-1);
	translator->addTranslation(164, 999, 1204, 3.51, 7.34e-4);
	translator->addTranslation(164, 1205, 1307, -6.37, 1.93e-1);
	translator->addTranslation(164, 1308, 1364, -6.81, 1.23e-1);
	//168
	translator->addTranslation(168, 600, 803, -9.58, 2.69e-1);
	translator->addTranslation(168, 804, 998, 22.5, 2.0e-1);
	translator->addTranslation(168, 999, 1204, 10.9, -1.38e-2);
	translator->addTranslation(168, 1205, 1307, -7.38, 3.64e-1);
	translator->addTranslation(168, 1308, 1364, -7.48, -3.79e-1);
	//172
	translator->addTranslation(172, 600, 803, -5.22, 6.23e-1);
	translator->addTranslation(172, 804, 998, 32.9, -1.534);
	translator->addTranslation(172, 999, 1204, -1.40e-1, -1.33);
	translator->addTranslation(172, 1205, 1307, -1.15, -1.22);
	translator->addTranslation(172, 1308, 1364, -1.09, -2.75e-1);
	config.addBuilder(translator);

	may2009DUTsOld1(config);
}

void BAToct2009(TbConfig &config) {
	cout << "Using configuration oct2009" << endl;
#ifdef SITEBATOCT2009_SET
	siteBATOct2009(config);
#else
	cout << "Site config for BAToct2009 not found" << endl;
#endif
	// Set up and configure eventbuilders
	config.addBuilder(new BatTrack);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	Chi2Builder* chi2 = new Chi2Builder(20.0); //Chi2 cut at 20.
	config.addBuilder((EventBuilder*) chi2);
	AngleCuts* angles = new AngleCuts("calibs/bat-oct2009-trackangle.txt", 1.5); //Angle cut at mean +/- 1.5 sigma
	config.addBuilder((EventBuilder*) angles);
	DutSync* dutsync = new DutSync();
	dutsync->readFile("calibs/Oct09_h8_Bat_all-checkdutsync-163-dutsync.txt");
	config.addBuilder(dutsync);

	DUT* m162 = new DUT("PLANAR", 162, 0);
	m162->lv1Range(3, 10);
	m162->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-162-masks.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus7degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus15degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus22degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus22degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus30degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus30degree-getetacorr-162-etacalib.txt");
	m162->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus7degree-getetacorr-162-etacalib.txt");
	config.addDut(m162);

	DUT* m163 = new DUT("STA", 163, 0);
	m163->lv1Range(3, 10);
	m163->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-163-masks.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus7degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus15degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus22degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus22degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus30degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus30degree-getetacorr-163-etacalib.txt");
	m163->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus7degree-getetacorr-163-etacalib.txt");
	config.addDut(m163);
	DUT* m165 = new DUT("FBK", 165, 0);
	m165->lv1Range(2, 10);
	m165->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-165-masks.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus7degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus15degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus22degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus22degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus30degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.minus30degree-getetacorr-165-etacalib.txt");
	m165->ecorrs.addEcorr(
			"calibs/Oct09_H8_Bat_FieldON.plus7degree-getetacorr-165-etacalib.txt");
	config.addDut(m165);
}

//Configuration for BAT Oct 2010
void BAToct2010(TbConfig &config) {
	cout << "Using configuration Oct 2010" << endl;
	//Site configuration
#ifdef SITEBATOCT2010_SET
	siteBATOct2010(config);
#else
	cout << "Site config for BAT  Oct 2010 not found" << endl;
#endif
	// Set up and configure eventbuilders
	config.addBuilder(new BatTrack);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	Chi2Builder* chi2 = new Chi2Builder(10.0); //Chi2 cut
	config.addBuilder((EventBuilder*) chi2);
	//AngleCuts* angles = new AngleCuts( "calibs/Oct2010-BAT-fieldOFF-0degree-trackangle.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
	AngleCuts* angles = new AngleCuts("calibs/Oct2010-BAT-trackangle.txt", 1.5); //Angle cut at mean +/- 1.5 sigma
	config.eventbuilders.push_back((EventBuilder*) angles);
	//DutSync* dutsync = new DutSync();
	//dutsync->readFile("calibs/dutsync-checkdutsync-160-dutsync.txt");
	//config.addBuilder( dutsync );

	//Configure DUTs used in the Oct 2009 test beam
	DUT* m160 = new DUT("STA-3G", 160, 0);
	//Range of acceptable lv1 values, keep loose to avoid false positives in noise scan
	m160->lv1Range(6, 10);
	m160->addMasks(
			"calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-160-masks.txt");
	//Eta corrections, output from getEtaCorr.
	//m160->addToTCalib("calibs/totToQ_160.root");
	//FIXME:
	//m160->totcalib->chargeCorrection(8.117/7.0); //wrong C-low
	config.addDut(m160);

	DUT* m161 = new DUT("FBK-3E-5E15", 161, 0);
	m161->lv1Range(6, 9);
	m161->addMasks(
			"calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-161-masks.txt");
	//m161->addToTCalib("calibs/totToQ_164.root");
	config.addDut(m161);

	DUT* m162 = new DUT("FBK-4E-3E15", 162, 0);
	m162->lv1Range(6, 10);
	m162->addMasks(
			"calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-162-masks.txt");
	//m162->addToTCalib("calibs/totToQ_168.root");
	config.addDut(m162);

	DUT* m163 = new DUT("PLANAR", 163, 0);
	m163->lv1Range(6, 9);
	m163->addMasks(
			"calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-163-masks.txt");
	//m163->addToTCalib("calibs/totToQ_172.root");
	config.addDut(m163);
}

void EudetPPS2009(TbConfig& config) {
#ifdef SITEEUDETOCT2009_SET
	siteEudetOct2009(config);
#else 
	cout << "Site config for Eudet2009 not found" << endl;
#endif
	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(14);
	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	DUT* m10 = new DUT("SIN", 10, 3);
	m10->lv1Range(0, 10);
	config.addDut(m10);
	DUT* m15 = new DUT("SIN", 15, 3);
	m15->lv1Range(0, 10);
	config.addDut(m15);
	DUT* m14 = new DUT("SIN", 14, 3);
	m14->lv1Range(0, 10);
	config.addDut(m14);
}

void SIMbat2010(TbConfig &config) {
	cout << "Using configuration SimBat2010" << endl;
#ifdef SITESIMBAT2010_SET
	siteSimBat2010(config);
#else
	cout << "Site config for simBat2010 not found" << endl;
#endif
	// Set up and configure eventbuilders
	/*config.eventbuilders.push_back(new BatTrack);
	 config.eventbuilders.push_back(new PixelMasker);
	 config.eventbuilders.push_back(new ClusterFinder());
	 Chi2Builder* chi2 = new Chi2Builder(15.0); //Chi2 cut at 20.
	 config.eventbuilders.push_back( (EventBuilder*) chi2);
	 //AngleCuts* angles = new AngleCuts( "calibs/bat-oct2009-trackangle.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
	 //config.eventbuilders.push_back( (EventBuilder*) angles );
	 //DutSync* dutsync = new DutSync();
	 //dutsync->readFile("calibs/Oct09_h8_Bat_all-checkdutsync-163-dutsync.txt");
	 //config.eventbuilders.push_back( dutsync );

	 config.eventbuilders.push_back(new simBaseBuilder());
	 //config.eventbuilders.push_back(new simTruthBuilder());
	 */

	//Data getters
	config.addBuilder(new BatTrack());
	config.addBuilder(new simBaseBuilder());
	config.addBuilder(new simTruthBuilder());
	config.addBuilder(new simPixelEdepBuilder());
	//Simulator
	config.addBuilder(new simDutRunner());
	//Masks & cuts
	config.addBuilder(new PixelMasker());
	config.addBuilder(new ClusterFinder());
	config.addBuilder(new CheckRegion());
	config.addBuilder(new Chi2Builder(15.0));
	//config.addBuilder(new etaCutter(0.2,0.8));

	DUT* m160 = new DUT("PLANAR", 160, 0);
	m160->lv1Range(0, 1);
	//m162->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-162-masks.txt");
	//m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-162-etacalib.txt");
	config.addDut(m160);
	config.simIdenNameMap[160] = "DUT1";

	DUT* m164 = new DUT("STA-3G", 164, 0);
	m164->lv1Range(0, 1);
	//m163->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-163-masks.txt");
	//m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-163-etacalib.txt");
	config.addDut(m164);
	config.simIdenNameMap[164] = "DUT2";

	DUT* m168 = new DUT("FBK-3E7", 168, 0);
	m168->lv1Range(0, 1);
	//m165->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-165-masks.txt");
	//m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-165-etacalib.txt");
	config.addDut(m168);
	config.simIdenNameMap[168] = "DUT3";

	DUT* m172 = new DUT("FBK-3EM5", 172, 0);
	m172->lv1Range(0, 1);
	//m165->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-165-masks.txt");
	//m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-165-etacalib.txt");
	config.addDut(m172);
	config.simIdenNameMap[172] = "DUT4";

	//Needed for simDataKeeper::getHitsByIden()
	config.simIdenNameMap_all[81] = "BatMod1";
	config.simIdenNameMap_all[83] = "BatMod3";
	config.simIdenNameMap_all[86] = "BatMod6";
	//Copy DUT idens into full list
	for (map<int, string>::iterator it = config.simIdenNameMap.begin();
			it != config.simIdenNameMap.end(); it++) {
		config.simIdenNameMap_all[it->first] = it->second;
	}
}

void rangeSupplied(bool gotRange, char* argc[]) {
	if (gotRange) {
		std::cerr << "More than one range of runs supplied" << std::endl;
		printUsage(argc);
	}
}

void listAnalysesAndExit(TbConfig &config){
	DUT *fakeDut = new DUT();
	config.analysisnames.clear();
	allAnalyses(config, fakeDut);
	cout << "List of available analyses:" << endl;
	for(list<TbAnalysis*>::const_iterator it=config.analysis.begin();
		  it!=config.analysis.end(); it++){
		cout << (*it)->basename << "\n";
		}
	cout << endl;
	exit(-1); //Must exit as the fake DUT would interfere with further analyses
}

void printUsage(char* argv[]) {
	cout
			<< "Needs ONE of the following arguments to determine which files to loop over:\n"
			<< " -l <runList> \t\t\t Loop over all runs in <runList>\n"
			<< " -r <firstrun> <lastrun> \t Loop over all runs from <firstrun> to <lastrun>\n"
			<< " -s <run> \t\t\t Loop over a single run\n\n"
			<< "Optional arguments:\n"
			<< " -a <analysisname> \t\t Only run analysis named <analysisname>. If not supplied, run over all available analyses (see -alist).\n"
			<< " -i <iden> \t\t\t Only run over iden <iden>. If not supplied, loop over all available idens \n"
			<< " -c <configname> \t\t Configures modules for configuration <configname>. If not supplied, default to that set in siteconfig.h\n"
			<< " -v <verbosity-level> \t\t Should be one of (in order of increasing verbosity) error, info, debug, debug2 or debug3. Default is info.\n"
			<< " -e <extension> \t\t Plot extension. Should be a file format recognized by ROOT.\n"
			<< " -d <data-path> \t\t Path to input tbtrackX.root files.\n"
			<< " -o <out-path>  \t\t Path to output directory.\n"
	                << " -org <organize-output>  \t Organizes output .root file into subdirectories by analysis name\n"
			<< " -b <base-name> \t\t Base name for all output. Defaults to name based on run list or range of runs.\n"
			<< " -f \t\t\t\t Redirect output(stdout) to log file named after the base name.\n"
			<< " -P:<key> <value> \t\t Pass an extra parameter, which can be used by analyses or builders.\n\n"
			<< " -alist \t\t\t Print a list of available analyses\n\n"
			<< "Simple Example:\n"
			<< " "
			<< argv[0]
			<< " -l runlists/may2009-bon-a15\n"
			<< "Applies all available analyses to all the runs in runlist for all idens.\n\n"
			<< "Complete Example:\n"
			<< " "
			<< argv[0]
			<< " -l runlists/may2009-bon-a15 -i 160 164 -a sumtot residuals -c may2009 -b test1 -e png -f -o out -P:param1 1.3 -P:param2 urk\n"
			<< "Applies analysis named \"sumtot\" and analysis named \"residuals\" to all runs in runlist for idens 160 and 164. "
			<< "Png plots are made and saved to relative path out. "
			<< "A log file(out/test1.log) is made, no output to screen. "
			<< "Extra-parameters \"param1\" and \"param2\" are stored internally with values \"1.3\" and \"urk\", and may be accessed from any builder and/or analysis."
			<< endl;
	exit(-1);
}

void parseArgs(int argc, char* argv[], TbConfig &config) {
	bool gotRange(false);
	int index = 0;
	while (index + 1 < argc) {
		index++;
		if (strcmp(argv[index], "-l") == 0 and argc > index + 1) { // runlist from file
			rangeSupplied(gotRange, argv);
			config.makeRunList(argv[++index]);
			gotRange = true;
			continue;
		}
		if (strcmp(argv[index], "-r") == 0 and argc > index + 2) { // range of runs
			int firstrun(atoi(argv[++index])), lastrun(atoi(argv[++index]));
			rangeSupplied(gotRange, argv);
			if (firstrun == 0 or lastrun == 0) {
				break;
			}
			config.makeRunList(firstrun, lastrun);
			gotRange = true;
			continue;
		}
		if (strcmp(argv[index], "-s") == 0 and argc > index + 1) { // single run
			int run = atoi(argv[++index]);
			rangeSupplied(gotRange, argv);
			if (run == 0) {
				break;
			}
			config.makeRunList(run);
			gotRange = true;
			continue;
		}
		if (strcmp(argv[index], "-v") == 0 and argc > index + 1) { // verbosity
			++index;
			if (strcmp(argv[index], "error") == 0) {
				config.logLevel = kERROR;
			} else if (strcmp(argv[index], "info") == 0) {
				config.logLevel = kINFO;
			} else if (strcmp(argv[index], "debug") == 0) {
				config.logLevel = kDEBUG;
			} else if (strcmp(argv[index], "debug2") == 0) {
				config.logLevel = kDEBUG2;
			} else if (strcmp(argv[index], "debug3") == 0) {
				config.logLevel = kDEBUG3;
			} else {
				cout << "Unknon verbosity level: " << argv[index] << endl;
				gotRange = false;
				break;
			}
			continue;
		}
		if (strcmp(argv[index], "-c") == 0 and argc > index + 1) { // config name
			config.tbslot = argv[++index];
			continue;
		}
		// -m option added to this file by Botho (probably copied from old version)
		if (strcmp(argv[index], "-m") == 0 and argc > index + 1) { // max number of runs
		  config.MaxNumOfTrigToBeProcPerRun = atoi(argv[++index]);
			continue;
		}
		if (strcmp(argv[index], "-o") == 0 and argc > index + 1) { //outPath
			sprintf(config.outPath, "%s", argv[++index]);
			continue;
		}
		if (strcmp(argv[index], "-d") == 0 and argc > index + 1) { //dataPath
			sprintf(config.dataPath, "%s", argv[++index]);
			continue;
		}
		if (strcmp(argv[index], "-e") == 0 and argc > index + 1) { //plotExtension
			sprintf(config.plotExtension, "%s", argv[++index]);
			continue;
		}
		if (strcmp(argv[index], "-b") == 0 and argc > index + 1) { //baseName
			sprintf(config.name, "%s", argv[++index]);
			continue;
		}
		if (strcmp(argv[index], "-i") == 0 and argc > index + 1) { // idens
			int iiden = index + 1;
			for (; iiden < argc; iiden++) {
				int iden = atoi(argv[iiden]);
				if (iden == 0) {
					break;
				}
				config.idens.push_back(iden);
			}
			index = iiden - 1;
			continue;
		}
		if (strcmp(argv[index], "-alist") == 0) { // list analyses
			listAnalysesAndExit(config);
			abort(); //listAnalysesAndExit has to quit tbmon
		}
		if (strcmp(argv[index], "-a") == 0 and argc > index + 1) { // analyses
			int iname = index + 1;
			for (; iname < argc; iname++) {
				char* name = (char*) argv[iname];
				if (name[0] == '-') {
					break;
				}
				config.analysisnames.push_back(name);
			}
			index = iname - 1;
			continue;
		}
		if (strcmp(argv[index], "-f") == 0) {
			config.logToFile = true;
			continue;
		}
		if (strcmp(argv[index], "-org") == 0) {
			config.organizeOutput = true;
			continue;
		}
		if (strncmp(argv[index], "-P:", 3) == 0 && argc > index + 1) {
			string key = argv[index];
			if (key.size() <= 3) {
				cerr
						<< "Error: Expected more after -P:, something like -P:<key> <val>."
						<< endl << endl;
				printUsage(argv);
			}
			key = key.substr(3);
			string val = argv[++index];
			//config.debug (maybe) not set yet
			//cout << "Found extraParam: " << key << " " << val << endl;

			//config.cmdLineExtras[key] = val;
			config.cmdLineExtras_tryset(key, val);

			continue;
		}
		std::cerr << "Malformed argument list." << endl;
		gotRange = false;
		break;
	}
	if (not gotRange) {
		printUsage(argv);
	}
	cout << "Preparing analysis " << config.name << endl;
}

int main(int argc, char* argv[]) {
	TbConfig config;
	//Parse args. See header
	parseArgs(argc, argv, config);
#ifdef SITECONFIG_SET
	siteConfig(config);
#else
	cout << "Global site config not found" << endl;
#endif
	// The configuration/calibration of DUT's and EventBuilders should be
	// performed in a function to make swiching configs as easy as possible.
	if (strcmp(config.tbslot, "may2009") == 0) {
		BATmay2009(config);
	} else if (strcmp(config.tbslot, "may2008") == 0) {
		BAT2008(config);
	} else if (strcmp(config.tbslot, "may2009old1") == 0) {
		BATmay2009old1(config);
	} else if (strcmp(config.tbslot, "batoct2009") == 0) {
		BAToct2009(config);
	} else if (strcmp(config.tbslot, "simBat2010") == 0) {
		SIMbat2010(config);
	} else if (strcmp(config.tbslot, "eudetoct2009") == 0) {
		EudetOct2009(config);
	} else if (strcmp(config.tbslot, "eudetnov2009") == 0) {
		EudetNov2009(config);
	} else if (strcmp(config.tbslot, "eudetpps2009") == 0) {
		EudetPPS2009(config);
	} else if (strcmp(config.tbslot, "eudetjune2010") == 0) {
		EudetJune2010(config);
	} else if (strcmp(config.tbslot, "eudetjuly2010") == 0) {
		EudetJuly2010(config);
	} else if (strcmp(config.tbslot, "batoct2010") == 0) {
		BAToct2010(config);
	} else if (strcmp(config.tbslot, "eudetnov2010") == 0) {
		EudetNov2010(config);
	} else if(strcmp(config.tbslot, "eudetfeb2011") == 0){
		EudetFeb2011(config);
	} else if (strcmp(config.tbslot, "eudetjune2011") == 0) {
		EudetJune2011(config);
	} else if (strcmp(config.tbslot, "eudetPPSjuly2011") == 0) {
		EudetPPSJuly2011(config);
	} else if (strcmp(config.tbslot, "eudetPPSSep2011") == 0) {
		EudetPPSSep2011(config);
	} else if (strcmp(config.tbslot, "CBjuly2013") == 0) {
		CBJuly2013(config);
	} else if (strcmp(config.tbslot, "eudetPPSMay2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSMay2012FEI4",config.name );
		EudetPPSMay2012FEI4(config);	
	} else if (strcmp(config.tbslot, "eudetPPSSep2012FEI3") == 0) {
		sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSSep2012FEI3",config.name );
		EudetPPSSep2012FEI3(config);
	} else if (strcmp(config.tbslot, "eudetPPSSep2012FEI4") == 0) {
		sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSSep2012FEI4",config.name );
		EudetPPSSep2012FEI4(config);	
	} else if (strcmp(config.tbslot, "eudetPPSMar2013FEI4") == 0) {
		sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSMar2013FEI4",config.name );
		EudetPPSMar2013FEI4(config);
	} else if (strcmp(config.tbslot, "eudetPPSAug2013FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSAug2013FEI4",config.name );
		EudetPPSAug2013FEI4(config);				
	} else if (strcmp(config.tbslot, "eudetPPSNov2013FEI34") == 0) {
		sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSNov2013FEI34",config.name );
		EudetPPSNov2013FEI34(config);
	} else if (strcmp(config.tbslot, "eudetPPSNov2013FEI4") == 0) {
		sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSNov2013FEI4",config.name );
		EudetPPSNov2013FEI4(config);
	} else if (strcmp(config.tbslot, "eudetPPSFeb2014FEI4") == 0) {
		sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSFeb2014FEI4",config.name );
		EudetPPSFeb2014FEI4(config);
	} else  if (strcmp(config.tbslot, "eudetITkOct2014FEI344") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetITkOct2014FEI344",config.name );
		EudetITkOct2014FEI344(config);
	} else  if (strcmp(config.tbslot, "eudetITKMar2015FEI4_25") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetITKMar2015FEI4_25",config.name );
		EudetITKMar2015FEI4_25(config);		
	} else  if (strcmp(config.tbslot, "eudetITKMar2015FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetITKMar2015FEI4",config.name );
		EudetITKMar2015FEI4(config);
	} else  if (strcmp(config.tbslot, "eudetIBLJun2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetIBLJun2012FEI4",config.name );
		EudetIBLJun2012FEI4(config);
	} else  if (strcmp(config.tbslot, "eudetIBLSep2011") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetIBLSep2011",config.name );
		EudetIBLSep2011(config);
	} else  if (strcmp(config.tbslot, "eudetAllpix") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetAllpix",config.name );
		EudetAllpix(config);
			
	} else {
		cout << "Test beam configuration called " << config.tbslot << " not found. Quitting" << endl;
		exit(-1);
	}

	//Graphic settings
	if (config.useAtlasStyle)
		SetAtlasStyle();
	else
		gROOT->SetStyle("Plain");
	gROOT->ProcessLine("gErrorIgnoreLevel = 20000;");
	gStyle->SetPalette(1);

	//Add the analysis objects
	for (list<DUT*>::iterator it = config.dutList.begin();
			it != config.dutList.end(); it++) {
		allAnalyses(config, (*it));
	}
	config.loop();
	return (0);
}
