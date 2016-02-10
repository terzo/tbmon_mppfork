#include "setupdriver.h"
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


//checks for existence of and creates setup directories if necessary
void setupDir(char* tosetup){
  struct stat st;
	if(stat(tosetup,&st)!=0){
	  char targetdir[800];
	  sprintf(targetdir,"mkdir -p %s", tosetup);
	  system(targetdir);
	  cout<<"New directory created."<<endl;
	}
}

void setupRunProc(TbConfig &config, void (*choiceSet)(TbConfig &, bool, bool, bool, char*) ){
  
  //Graphic settings
  if (config.useAtlasStyle)
    SetAtlasStyle();
  else
    gROOT->SetStyle("Plain");
  gROOT->ProcessLine("gErrorIgnoreLevel = 20000;");
  gStyle->SetPalette(1);
  
  char outdir[800];
  strcpy(outdir, (const char*) config.outPath);
  
  //hotpixelfinder setup
  sprintf(config.outPath, "%s%s/", outdir,(char*) "hotpixelfinder" );
  setupDir(config.outPath);
  choiceSet(config, false, false, false, (char*) outdir);
  for (list<DUT*>::iterator it = config.dutList.begin();
       it != config.dutList.end(); it++) {
    config.addAnalysis(new HotPixelFinder, "hotpixelfinder", (*it));
  }
  config.loop();
  config.clear();
  
  //first checkalign setup
  sprintf(config.outPath, "%s%s/", outdir,(char*) "checkalign1" );
  setupDir(config.outPath);
  choiceSet(config, true, false, false, (char*) outdir);
  for (list<DUT*>::iterator bit = config.dutList.begin();
       bit != config.dutList.end(); bit++) {
    config.addAnalysis(new CheckAlign, "checkalign", (*bit));
  }
  config.loop();
  config.clear();
  
  //getetacorr setup
  sprintf(config.outPath, "%s%s/", outdir,(char*) "getetacorr" );
  setupDir(config.outPath);
  choiceSet(config, true, true, false, (char*) outdir);
  for (list<DUT*>::iterator it = config.dutList.begin();
       it != config.dutList.end(); it++) {
    config.addAnalysis(new GetEtaCorr, "getetacorr", (*it));
  }
  config.loop();
  config.clear();
  
  //second checkalign setup
  sprintf(config.outPath, "%s%s/", outdir,(char*) "checkalign2" );
  setupDir(config.outPath);
  choiceSet(config, true, false, true, (char*) outdir);
  for (list<DUT*>::iterator it = config.dutList.begin();
       it != config.dutList.end(); it++) {
    config.addAnalysis(new CheckAlign, "checkalign", (*it));
  }
  config.loop();

}

void EudetIBLSep2011(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
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
	config.addBuilder(new Chi2Builder(20.0));
	
	
	// mask noisy pixels?
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;

	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;	

	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	dutnumber = 22;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 23;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 500.0, 25.0, 10.0,10.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0, 5.0, 5.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);


dutnumber = 22;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m22 = makeFEI4("SCM132", dutnumber, 2, 500.0, 25.0, 10.0,10.0);
	m22->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m22->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m22->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m22->addToTCalib(calibfile.c_str());
	}
	config.addDut(m22);

	dutnumber = 23;
	DUT* m23 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0, 5.0, 5.0);
	m23->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m23->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m23->addMasks(maskfile);
	}
	if (addtotcalib) {
		m23->addToTCalib(calibfile.c_str());
	}
	config.addDut(m23);
}



void EudetAllpix(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;	

	
	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 500.0, 25.0, 10.0,10.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0, 5.0, 5.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);
/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}

void EudetIBLJun2012FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETIBLJUN2012FEI4_SET
	siteEudetIBLJun2012FEI4(config);
#else 
	cout << "Site config for EudetITKMar2015FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;	

	
	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	
	dutnumber = 22;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 500.0, 25.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);
	
		dutnumber = 22;
	DUT* m22 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m22->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m22->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m22->addMasks(maskfile);
	}
	if (addtotcalib) {
		m22->addToTCalib(calibfile.c_str());
	}
	config.addDut(m22);
/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}


void EudetITKMar2015FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;	

	
	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 500.0, 25.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);
/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}

void EudetITKMar2015FEI4_25(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;	

	
	dutnumber = 20;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);

	dutnumber = 21;
	sprintf(alignfile, "%scheckalign%d/%s-checkalign-%d-translations.txt",
		config.outPath, calibdir, config.name, dutnumber);
	cout << "Applying translations " << alignfile << " for DUT "
			<< dutnumber << endl;
	eut->addTranslation(alignfile, dutnumber);
	

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4_25("SCM132", dutnumber, 2, 500.0, 25.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);
/*
	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
*/
}


void EudetPPSFeb2014FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSFEB2014FEI4_SET
	siteEudetPPSFeb2014FEI4(config);
#else 
	cout << "Site config for eudetPPSFeb2014FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


void EudetPPSNov2013FEI34(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSNOV2013FEI34_SET
	siteEudetPPSNov2013FEI34(config);
#else 
	cout << "Site config for eudetPPSNov2013FEI34 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("VTT", dutnumber, 2, 600.0, 50.0);
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m10->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 20;
	DUT* m20 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


void EudetPPSAug2013FEI34(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSAUG2013FEI34_SET
	siteEudetPPSAug2013FEI34(config);
#else 
	cout << "Site config for eudetPPSAug2013FEI34 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("VTT", dutnumber, 2, 600.0, 50.0);
	//DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m10->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 20;
	DUT* m20 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSAug2013FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSAUG2013FEI4_SET
	siteEudetPPSAug2013FEI4(config);
#else 
	cout << "Site config for eudetPPSAug2013 FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSMar2013FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSMAR2013FEI4_SET
	siteEudetPPSMar2013FEI4(config);
#else 
	cout << "Site config for eudetPPSMar2013 FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


void EudetPPSSep2012FEI3(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("VTT", dutnumber, 2, 600.0, 50.0);
	m10->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);
	
	dutnumber = 11;
	DUT* m11 = makeFEI3("VTT", dutnumber, 2, 600.0, 50.0);
	m11->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m11->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m11->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m11->addToTCalib(calibfile.c_str());
	}
	config.addDut(m11);
	
	dutnumber = 12;
	DUT* m12 = makeFEI4("SCM131", dutnumber, 2, 600.0, 50.0);
	m12->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m12->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m12->addMasks(maskfile);
	}
	if (addtotcalib) {
		m12->addToTCalib(calibfile.c_str());
	}
	config.addDut(m12);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 11;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);
		
		dutnumber = 12;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;		
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}
/*
void EudetIBLJun2012FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETIBLJUN2012FEI4_SET
	siteEudetIBLJun2012FEI4(config);
#else 
	cout << "Site config for eudetIBLJun2012 FEI4 not found" << endl;
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
	config.addBuilder(new Chi2Builder(100.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("VTT10", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("CNM55", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}
*/
void EudetIBLJun2012FEI4_3DUT(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETIBLJUN2012FEI4_3DUT_SET
	siteEudetIBLJun2012FEI4_3DUT(config);
#else 
	cout << "Site config for eudetIBLJun2012_3DUT FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("FBK13", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 14); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("FBK11", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 14); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);
	
	dutnumber = 22;
	DUT* m22 = makeFEI4("CNM55", dutnumber, 2, 250.0, 50.0);
	m22->lv1Range(0, 14); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m22->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m22->addMasks(maskfile);
	}
	if (addtotcalib) {
		m22->addToTCalib(calibfile.c_str());
	}
	config.addDut(m22);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 22;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);
		
		config.addBuilder(translator);
	}
}

void EudetPPSSep2012FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSSEP2012FEI4_SET
	siteEudetPPSSep2012FEI4(config);
#else 
	cout << "Site config for eudetPPSSep2012 FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM132", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM131", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSMay2012FEI4(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSMAY2012FEI4_SET
	siteEudetPPSMay2012FEI4(config);
#else 
	cout << "Site config for eudetPPSMay2012 FEI4 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 20;
	//DUT* m10 = makeFEI3("DO-I-7", dutnumber, 2, 600.0, 50.0);
	DUT* m20 = makeFEI4("SCM127", dutnumber, 2, 250.0, 50.0);
	m20->lv1Range(3, 9); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m20->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m20->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m20->addToTCalib(calibfile.c_str());
	}
	config.addDut(m20);

	dutnumber = 21;
	DUT* m21 = makeFEI4("SCM130", dutnumber, 2, 250.0, 50.0);
	m21->lv1Range(3, 9); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m21->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m21->addMasks(maskfile);
	}
	if (addtotcalib) {
		m21->addToTCalib(calibfile.c_str());
	}
	config.addDut(m21);

	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 20;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 21;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSMay2012(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSMAY2012_SET
	siteEudetPPSMay2012(config);
#else 
	cout << "Site config for eudetPPSMay2012 not found" << endl;
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
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("LUB1e16", dutnumber, 2, 600.0, 50.0);
	//DUT* m20 = makeFEI4("SCM127", dutnumber, 2, 250.0, 50.0);
	m10->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 11;
	DUT* m11 = makeFEI3("SLID10", dutnumber, 2, 600.0, 50.0);
	m11->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m11->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m11->addMasks(maskfile);
	}
	if (addtotcalib) {
		m11->addToTCalib(calibfile.c_str());
	}
	config.addDut(m11);

	dutnumber = 12;
	DUT* m12 = makeFEI3("SLID3", dutnumber, 2, 600.0, 50.0);
	m12->lv1Range(0, 16); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m12->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m12->addMasks(maskfile);
	}
	if (addtotcalib) {
		m12->addToTCalib(calibfile.c_str());
	}
	config.addDut(m12);



	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 11;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 12;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}

void EudetPPSMay2012_2duts(TbConfig &config, bool hotpix, bool aptrans, bool etacor, char* calibdir) {
#ifdef SITEEUDETPPSMAY2012_2duts_SET
	siteEudetPPSMay2012_2duts(config);
#else 
	cout << "Site config for eudetPPSMay2012 not found" << endl;
#endif

	sprintf(config.treeName, "eutracks");
	EuBuildTrack* eut = new EuBuildTrack();
	eut->addMatchDUT(10);
	eut->addMatchDUT(11);

	eut->nMatches(1);
	config.addBuilder(eut);
	config.addBuilder(new PixelMasker);
	config.addBuilder(new ClusterFinder);
	config.addBuilder(new CheckRegion);
	config.addBuilder(new Chi2Builder(15.0));

	// mask noisy pixels?
	bool masknoisyanddeadpixels = hotpix;
	// apply eta corrections?
	bool applyetacorr = etacor;
	// add TOT vs charge calibration?
	bool addtotcalib = false;
	// apply translation to center residuals at zero
	bool applytranslation = aptrans;


	char etacorrfile[300];
	char maskfile[300];
	char alignfile[300];
	std::string calibfile = "/cal.out";

	int dutnumber;

	dutnumber = 10;
	DUT* m10 = makeFEI3("LUB1e16", dutnumber, 2, 600.0, 50.0);
	//DUT* m20 = makeFEI4("SCM127", dutnumber, 2, 250.0, 50.0);
	m10->lv1Range(3, 9); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m10->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m10->addMasks(maskfile);
	}
	if (addtotcalib) {
		//Add the analysis objects
		m10->addToTCalib(calibfile.c_str());
	}
	config.addDut(m10);

	dutnumber = 11;
	DUT* m11 = makeFEI3("SLID3", dutnumber, 2, 600.0, 50.0);
	m11->lv1Range(3, 9); // these cuts are only needed to look for hot pixels (in out-of-time events)
	if (applyetacorr) {
		// adding eta correction file. You may get it running the analysis getetacorr
		sprintf(etacorrfile, "%sgetetacorr/%s-getetacorr-%d-etacalib.txt",
				calibdir, config.name, dutnumber);
		cout << "Applying eta corrections " << etacorrfile << " for DUT "
				<< dutnumber << endl;
		m11->ecorrs.addEcorr(etacorrfile);
	}
	if (masknoisyanddeadpixels) {
		// adding noisy/dead pixel file. You may get it running the analysis hotpixelfinder
		sprintf(maskfile, "%shotpixelfinder/%s-hotpixelfinder-%d-masks.txt", calibdir,
				config.name, dutnumber);
		cout << "Applying hot pixel mask " << maskfile << " for DUT "
				<< dutnumber << endl;
		m11->addMasks(maskfile);
	}
	if (addtotcalib) {
		m11->addToTCalib(calibfile.c_str());
	}
	config.addDut(m11);



	Translator* translator = 0;
	if (applytranslation) {

		translator = new Translator();

		dutnumber = 10;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		dutnumber = 11;
		sprintf(alignfile, "%scheckalign1/%s-checkalign-%d-translations.txt",
			calibdir, config.name, dutnumber);
		cout << "Applying translations " << alignfile << " for DUT "
				<< dutnumber << endl;
		translator->addTranslation(alignfile, dutnumber);

		config.addBuilder(translator);
	}
}


void rangeSupplied(bool gotRange, char* argc[]) {
	if (gotRange) {
		std::cerr << "More than one range of runs supplied" << std::endl;
		printUsage(argc);
	}
}

void printUsage(char* argv[]) {
	cout
			<< "Needs ONE of the following arguments to determine which files to loop over:"
			<< endl
			<< " -l <runList> \t\t\t Loop over all runs in <runList>"
			<< endl
			<< " -r <firstrun> <lastrun> \t Loop over all runs from <firstrun> to <lastrun>"
			<< endl
			<< " -s <run> \t\t\t Loop over a single run"
			<< endl
			<< endl
			<< "Optional arguments:"
			<< endl
			<< " -a <analysisname> \t\t Only run analysis named <analysisname>. If not supplied, run over all available analyses."
			<< endl
			<< " -i <iden> \t\t\t Only run over iden <iden>. If not supplied, loop over all available idens."
			<< endl
			<< " -c <configname> \t\t Configures modules for configuration <configname>. If not supplied, default to that set in siteconfig.h"
			<< endl
			<< " -v <verbosity-level> \t\t Should be one of (in order of increasing verbosity) error, info, debug, debug2 or debug3. Default is info."
			<< endl
			<< " -e <extension> \t\t Plot extension. Should be a file format recognized by ROOT."
			<< endl
			<< " -d <data-path> \t\t Path to input tbtrackX.root files."
			<< endl
			<< " -o <out-path>  \t\t Path to output directory."
			<< endl
			<< " -b <base-name> \t\t Base name for all output. Defaults to name based on run list or range of runs."
			<< endl
			<< " -f \t\t\t\t Redirect output(stdout) to log file named after the base name."
			<< endl
			<< " -P:<key> <value> \t\t Pass an extra parameter, which can be used by analyses or builders."
			<< endl
			<< endl
			<< "Simple Example:"
			<< endl
			<< " "
			<< argv[0]
			<< " -l runlists/may2009-bon-a15"
			<< endl
			<< "Applies all available analyses to all the runs in runlist for all idens."
			<< endl
			<< endl
			<< "Complete Example:"
			<< endl
			<< " "
			<< argv[0]
			<< " -l runlists/may2009-bon-a15 -i 160 164 -a sumtot residuals -c may2009 -b test1 -e png -f -o out -P:param1 1.3 -P:param2 urk"
			<< endl
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
	if (strcmp(config.tbslot, "eudetPPSmay2012") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSmay2012",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSMay2012);
	} else if (strcmp(config.tbslot, "eudetPPSmay2012_2duts") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSmay2012_2duts",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSMay2012_2duts);
	} else if (strcmp(config.tbslot, "eudetPPSMay2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSmay2012FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSMay2012FEI4);
	} else if (strcmp(config.tbslot, "eudetPPSSep2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSSep2012FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSSep2012FEI4);
	} else if (strcmp(config.tbslot, "eudetIBLJun2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetIBLJun2012FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetIBLJun2012FEI4);
	} else if (strcmp(config.tbslot, "eudetPPSSep2012FEI3") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSSep2012FEI3",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSSep2012FEI3);
	} else if (strcmp(config.tbslot, "eudetPPSSep2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSSep2012FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSSep2012FEI4); 
	} else if (strcmp(config.tbslot, "eudetIBLJun2012FEI4_3DUT") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetIBLJun2012FEI4_3DUT",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetIBLJun2012FEI4_3DUT);	
	} else if (strcmp(config.tbslot, "eudetPPSMar2013FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSMar2013FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSMar2013FEI4);
	} else if (strcmp(config.tbslot, "eudetPPSAug2013FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSAug2013FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSAug2013FEI4);
	} else if (strcmp(config.tbslot, "eudetPPSAug2013FEI34") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSAug2013FEI34",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSAug2013FEI34);	
	} else if (strcmp(config.tbslot, "eudetPPSNov2013FEI34") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSNov2013FEI34",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSNov2013FEI34);
	} else if (strcmp(config.tbslot, "eudetPPSFeb2014FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetPPSFeb2014FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetPPSFeb2014FEI4);
	} else if (strcmp(config.tbslot, "eudetITKMar2015FEI4_25") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetITKMar2015FEI4_25",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetITKMar2015FEI4_25); 
	} else if (strcmp(config.tbslot, "eudetITKMar2015FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetITKMar2015FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetITKMar2015FEI4); 
	} else if (strcmp(config.tbslot, "eudetIBLJun2012FEI4") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetIBLJun2012FEI4",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetIBLJun2012FEI4); 	
	} else if (strcmp(config.tbslot, "eudetAllpix") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetAllpix",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetAllpix); 	
	} else if (strcmp(config.tbslot, "eudetIBLSep2011") == 0) {
	  sprintf(config.outPath, "%s%s/%s/", config.outPath, (char*) "eudetIBLSep2011",config.name );
	  setupDir(config.outPath);
	  setupRunProc(config, EudetIBLSep2011); 
	} else {
		cout << "Test beam configuration called " << config.tbslot << " not found. Quitting" << endl;
		exit(-1);
	}
	return (0);

}
