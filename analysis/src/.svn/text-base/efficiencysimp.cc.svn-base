#include "efficiencysimp.h"

void EfficiencySimp::readMask(string maskpath, TbConfig & config) {
	maskfile.open(maskpath.c_str());
	if (!maskfile.is_open()) {
		TBALOG(kINFO) << "Maskfile is not open" << endl;
		maskon = false;
		return;
	}
	maskfile.seekg(0, ios_base::beg);
	char* readin = new char[18];
	for (int r = 0; r < 160; r++) {
		maskfile.getline(readin, 19);
		for (int c = 0; c < 18; c++) {
			TBALOG (kDEBUG) << "Reading in mask" << endl;
			TBALOG (kDEBUG) << c << " " << readin[c] << " ";
			if (readin[c] == '1')
				mask[r][c] = true;
			else if (readin[c] == '0')
				mask[r][c] = false;
			else {
				TBALOG(kINFO) << "Maskfile badly formatted " << r << " " << c
						<< " " << readin[c] << endl;
				maskon = false;
				return;
			}
		}
		while (maskfile.peek() < 33 && !maskfile.eof()) {
			maskfile.ignore(1);
		}
	}
}

void EfficiencySimp::init(TbConfig & config) {

	//read in mask file, if present
	/*
	 string maskpath = config.cmdLineExtras_argGetter("A_efficiencysimp_maskfile",string(""),"Mask file path");
	 if (maskpath=="") maskon = false;
	 else {
	 maskon = true;
	 mask = new bool* [160];
	 for (int r = 0;r < 160; r++) {mask[r] = new bool [18];}
	 readMask(maskpath,config);
	 }
	 */
	maskon = false;

	//lvl1cut parameters
	/*
	 string lvl1string = config.cmdLineExtras_argGetter("A_efficiencysimp_lvl1values",string(""),"Lvl1 values");
	 if (lvl1string=="") lvl1cut = false;
	 else {
	 lvl1cut = true;
	 string ints;
	 stringstream lvl1stream(lvl1string);
	 while(getline(lvl1stream,ints,' ')) {
	 lvl1.push_back(atoi(ints.c_str()));
	 }
	 TBALOG (kDEBUG) << "Reading in lvl1 cut values: ";
	 for(vector<int>::const_iterator i = lvl1.begin(); i != lvl1.end(); i++) {TBALOG (kDEBUG) << (*i) << " ";}
	 }
	 */
	lvl1cut = false;

	//check for charge conversion parameters

	//list<DUT*>::const_iterator i;
	//for(i = config.dutList.begin(); i != config.dutList.end(); i++) {if (iden == (*i)->iden) break;}
	//if((*i)->totcalib != NULL) chargeconversion = true;
	//else chargeconversion = false;
	chargeconversion = false;

	/*
	 for (int c = 0; c < 18; c++) {
	 for (int r = 0; r < 160; r++) {
	 cout << calC[c][r] << " ";
	 }
	 }
	 */

	ntracks = 0;

	hitMap = new TH2D("", "Hit Map;Column;Row", 18, 0, 18, 160, 0, 160);
	anyhitMap = new TH2D("", "Any Hit Map;Column;Row", 18, 0, 18, 160, 0, 160);

	totHist = new TH1I("", "TOT dist", 256, -0.5, 255.5);
	totMap = new TH2D("", "Mean TOT Map", 18, 0, 18, 160, 0, 160);
	chargeClusterTree = new TTree("chargeclustertree", "chargeclustertree");
	chargeClusterTree->Branch("charge", &chargesum);
	lvl1Hist = new TH1I("", "Lvl1 dist", 16, 0, 16);
	chargeClusterMap = new TH2D("", "Mean cluster charge map", 18, 0, 18, 160,
			0, 160);

	//initialise tot collections
	for (int i = 0; i < 18; i++) {
		for (int j = 0; j < 160; j++) {
			totStore[i][j] = new vector<int>;
		}
	}
	for (int i = 0; i < 18; i++) {
		for (int j = 0; j < 160; j++) {
			chargeClusterStore[i][j] = new vector<float>;
		}
	}

}

void EfficiencySimp::event(const TbConfig &config, const Event &event) {

	ntracks++;

	TBALOG( kDEBUG2 ) << "Event: " << config.currentEntry << " Run: "
			<< config.currentRun
			<< endl;
	TBALOG( kDEBUG2 ) << "Accepted track." << endl;

	bool cut;
	//plot the position of *any* single hits and store TOT values
	for (vector<PllHit*>::const_iterator i = event.hits.begin();
			i != event.hits.end(); ++i) {
		if ((*i)->iden != this->iden)
			continue;
		if (maskon && !mask[(*i)->row][(*i)->col])
			continue;
		cut = true;
		if (lvl1cut) {
			for (vector<int>::const_iterator it = lvl1.begin();
					it != lvl1.end(); it++) {
				if ((*it) == (*i)->lv1)
					cut = false;
			}
		}
		if (lvl1cut && cut)
			continue;
		anyhitMap->Fill((*i)->col, (*i)->row);
		totHist->AddBinContent((*i)->tot);
		(totStore[(*i)->col][(*i)->row])->push_back((*i)->tot);
	}

	int coli, rowi, toti;
	//plot charge seen for each cluster
	if (chargeconversion) {
		for (vector<event::cluster_t>::const_iterator i =
				event.clusters.begin(); i != event.clusters.end(); ++i) {

			toti = 0;
			for (vector<PllHit*>::const_iterator j = (i)->begin();
					j != (i)->end(); ++j) {
				if ((*j)->iden != this->iden)
					continue;
				if ((*j)->tot > toti) {
					toti = (*j)->tot;
					coli = (*j)->col;
					rowi = (*j)->row;
				}
			}
			chargesum = cluster::getSumChargePP((*i), event);
			chargeClusterTree->Fill();
			if (coli > 17 || rowi > 159)
				continue;
			(chargeClusterStore[coli][rowi])->push_back(chargesum);
		}
	}

}

void EfficiencySimp::finalize(const TbConfig &config) {
//   TStyle* style = (TStyle*)gStyle->Clone("curstyle");
//   tbutils::atlasHistStyle();
	totHist->SetObjectStat(true);

	config.saveToFile(name, "hitsany", anyhitMap);

	config.saveToFile(name, "totdist", totHist);

	if (chargeconversion) {
		//TH1I *htemp = new TH1I("charge","charge",10000,0.0,10000.0);
		chargeClusterTree->Draw("charge");
		TH1I* htemp = (TH1I*) gPad->GetPrimitive("htemp");
		config.saveToFile(name, "chargeclusterdist", htemp);
	}
	TBALOG( kINFO ) << " Tracks:             " << ntracks << endl;
	TBALOG( kINFO ) << " Maximum tot bin:    "
			<< (totHist->GetMaximumBin() - 1) << endl;

	//compute ToT average for each pixel and enter in histogram
	int nentv;
	double aver;
	for (int i = 0; i < 18; i++) {
		for (int j = 0; j < 160; j++) {
			aver = 0.0;
			nentv = (totStore[i][j])->size();
			for (int it = 0; it < nentv; it++) {
				aver += (totStore[i][j])->at(it);
			}
			aver /= nentv;
			totMap->SetBinContent(i + 1, j + 1, aver);
		}
	}

	config.saveToFile(name, "ToTmap", totMap);

	/*vector<PllHit*> cut;
	 cut.reserve(lvl1Hist->GetBinContent(floor(cutm)) + lvl1Hist->GetBinContent(ceil(cutm)));
	 for (vector<PllHit*>::const_iterator i = event.hits.begin();i!=event.hits.end();i++) {
	 if ((*i)->lv1 == floor(cutm) || (*i)->lv1 == ceil(cutm)) cut.push_back(*i);
	 }
	 
	 event.hits = cut;*/

	if (chargeconversion) {
		for (int i = 0; i < 18; i++) {
			for (int j = 0; j < 160; j++) {
				aver = 0.0;
				nentv = (chargeClusterStore[i][j])->size();
				cout << nentv << endl;
				for (int it = 0; it < nentv; it++) {
					aver += (chargeClusterStore[i][j])->at(it);
				}
				aver /= nentv;
				chargeClusterMap->SetBinContent(i + 1, j + 1, aver);
			}
		}
	}

	if (chargeconversion)
		config.saveToFile(name, "chargeclustermap", chargeClusterMap);

}

