#include "eubuildtrack.h"
#include "dut.h"

void EuBuildTrack::init(TbConfig& config) {
	//Track stuff
	t_posX = new vector<double>();
	t_posY = new vector<double>();
	t_dxdz = new vector<double>();
	t_dydz = new vector<double>();
	t_chi2 = new vector<double>();
	t_iden = new vector<int>();
	t_trackNum = new vector<int>();
	t_ndof = new vector<double>();

	//Pixel stuff
	p_col = new vector<int>();
	p_row = new vector<int>();
	p_tot = new vector<int>();
	p_iden = new vector<int>();
	p_lv1 = new vector<int>();
	p_chip = new vector<int>();

	//global rotation tree
	gr_ID = new std::vector<int>;
	gr_Alpha = new std::vector<double>;
	gr_Beta  = new std::vector<double>;
	gr_Gamma = new std::vector<double>;
	gr_RotXY = new std::vector<double>;
	gr_RotZX = new std::vector<double>;
	gr_RotZY = new std::vector<double>;
	gr_RotXYErr = new std::vector<double>;
	gr_RotZXErr = new std::vector<double>;
	gr_RotZYErr = new std::vector<double>;


	for (list<int>::iterator it = m_idens.begin(); it != m_idens.end(); it++) {
		m_hitMap[(*it)] = vector<PllHit*>();
	}
}

/**
 * set addresses of tbtrack file ttree to those of internal eubuildtrack ttree pointers
 * at the moment only considers eutracks and zspix branches
 */
void EuBuildTrack::initRun(TbConfig& config) {
	m_translations.clear();
	for (list<TransCalib>::iterator it = m_calibs.begin(); it != m_calibs.end();
			it++) {
		if (config.currentRun >= (*it).firstRun
				and config.currentRun <= (*it).lastRun) {
			m_translations[(*it).iden] = (*it);
		}
	}
	tracktree = (TTree*) config.infile->Get("eutracks");

	tracktree->SetBranchAddress("nTrackParams", &t_nTrackParams);
	tracktree->SetBranchAddress("euEvt", &t_euEv);
	tracktree->SetBranchAddress("xPos", &t_posX);
	tracktree->SetBranchAddress("yPos", &t_posY);
	tracktree->SetBranchAddress("dxdz", &t_dxdz);
	tracktree->SetBranchAddress("dydz", &t_dydz);
	tracktree->SetBranchAddress("iden", &t_iden);
	tracktree->SetBranchAddress("trackNum", &t_trackNum);
	tracktree->SetBranchAddress("chi2", &t_chi2);
	tracktree->SetBranchAddress("ndof", &t_ndof);

	pixeltree = (TTree*) config.infile->Get("zspix");
	pixeltree->SetBranchAddress("nPixHits", &p_nHits);
	pixeltree->SetBranchAddress("col", &p_col);
	pixeltree->SetBranchAddress("row", &p_row);
	pixeltree->SetBranchAddress("tot", &p_tot);
	pixeltree->SetBranchAddress("lv1", &p_lv1);
	pixeltree->SetBranchAddress("chip", &p_chip);
	pixeltree->SetBranchAddress("iden", &p_iden);
	pixeltree->SetBranchAddress("euEvt", &p_euEv);
	p_euEv = 0;

	//config.infile->GetObject("DUTrotation", globalrotationtree);
	globalrotationtree = (TTree*) config.infile->Get("DUTrotation");
	if (globalrotationtree) {
		hasRotationTree = true;
		// rotation information available - will be read in now
		cout << "Data file " << config.currentRun
				<< " contains rotation tree. It will be now read." << endl;
		globalrotationtree->SetBranchAddress("ID", &gr_ID);
		globalrotationtree->SetBranchAddress("Alpha", &gr_Alpha);
		globalrotationtree->SetBranchAddress("Beta", &gr_Beta);
		globalrotationtree->SetBranchAddress("Gamma", &gr_Gamma);
		globalrotationtree->SetBranchAddress("RotXY", &gr_RotXY);
		globalrotationtree->SetBranchAddress("RotZX", &gr_RotZX);
		globalrotationtree->SetBranchAddress("RotZY", &gr_RotZY);
		globalrotationtree->SetBranchAddress("RotXYErr", &gr_RotXYErr);
		globalrotationtree->SetBranchAddress("RotZXErr", &gr_RotZXErr);
		globalrotationtree->SetBranchAddress("RotZYErr", &gr_RotZYErr);
		
		int treeEntries = globalrotationtree->GetEntries();
		int idSize = gr_ID->size();

		if (treeEntries < 1) {
			TBALOG(kERROR) << "No angle collections found" << endl;
		} else if (gr_Alpha->size() != idSize || gr_RotXY->size() != idSize
				|| gr_Beta->size() != idSize || gr_RotZX->size() != idSize
				|| gr_Gamma->size() != idSize || gr_RotZY->size() != idSize) {
			TBALOG(kERROR) << "Corrupted rotation data. Vector size mismatch."
					<< endl;
		} else {
			if (treeEntries > 1) {
				TBALOG(kINFO) << "More than one angle collection found. "
						<< "Using first only." << endl;
			}
			globalrotationtree->GetEntry(0);
			for (int id = 0; id < idSize; id++) {
				//Check if DUT exists and get pointer
				DUT* dut = 0;
				map<int, DUT*>::const_iterator it = config.dutMap.find(
						gr_ID->at(id));
				if (it == config.dutMap.end()) {
					continue;
				} else {
					dut = (*it).second;
				}

				double degToRad = TMath::Pi() / 180;
				dut->anglePhiFromGEAR = degToRad * gr_RotZY->at(id);
				dut->angleEtaFromGEAR = degToRad * gr_RotZX->at(id);
				dut->anglePhiFromAlignment = degToRad * gr_Alpha->at(id);
				dut->angleEtaFromAlignment = degToRad * gr_Beta->at(id);
				dut->hasAnglesFromReco = true;
			}
		}
	} else {
		hasRotationTree = false;
		// rotation tree not available in this data file
		cout << "Data file " << config.currentRun
				<< " does not contain rotation tree." << endl;
	}
}

void EuBuildTrack::initEvent(TbConfig &config, map<int, Event> &events) {
	tracktree->GetEntry(config.currentEntry);
	pixeltree->GetEntry(config.currentEntry);
	if (getHasRotationTree() == true) {
		globalrotationtree->GetEntry(0);
	};
	//sort hits
	pllCounter = 0;
	for (map<int, std::vector<PllHit*> >::iterator it = m_hitMap.begin();
			it != m_hitMap.end(); it++) {
		(*it).second.clear();
	}
	for (int hit = 0; hit < p_iden->size(); hit++) {
		//We only want hits for events and match DUTs
		bool inEvent = (events.find(p_iden->at(hit)) != events.end());
		bool matchDut = (m_hitMap.find(p_iden->at(hit)) != m_hitMap.end());
		if ((not inEvent) and (not matchDut)) {
			continue;
		}
		//Waste less memory, rewrite instead of reallocate
		PllHit* apix;
		if (pllCounter < allPllHits.size()) {
			apix = allPllHits.at(pllCounter);
		} else {
			apix = new PllHit();
			allPllHits.push_back(apix);
		}
		pllCounter++;

		apix->iden = p_iden->at(hit);
		apix->col = p_col->at(hit);
		apix->row = p_row->at(hit);
		/*if(apix->iden==20)
		{  
			if (apix->col % 2!=0)
  			{
  			  apix->row *= 2;
			  apix->row++;
  			  apix->col = (apix->col-1)/2;
  			}
  			else 
  			{
  			  apix->row *= 2;
  			  apix->col = apix->col/2;
  			}
		}*/
		apix->tot = p_tot->at(hit);
		apix->lv1 = p_lv1->at(hit);
		apix->trig = p_euEv;
		apix->chp = p_chip->at(hit);
		if (matchDut) {
			m_hitMap[apix->iden].push_back(apix);
		}
		if (inEvent) {
			events[apix->iden].rawHits.push_back(apix);
		}
	}
}

/**
 * Loop over assigned matchDUTs and look for matching tracks in the ones, that
 * the current event does not belong to.
 * Then assign the track parameters of the matched track at the position of the
 * DUT of the event to this event.
 */
void EuBuildTrack::buildEvent(Event &event, map<int, Event> &,
		TbConfig & config) {
	//event.fBase = event::kGood;
	//event.fTrack = event::kBad;
	double transX(0.0), transY(0.0);
	double trackNum = -1;
	// double chi2Min = 999; changed by Botho
	double chi2Min = 999999;
	event.track = NULL;
	bool first = true;
	int nMatch = 0;
	//Loop over all track idens
	for (list<int>::iterator it = m_idens.begin(); it != m_idens.end(); it++) {
		//Only look at "other" DUTs, not the one which Event belongs to.
		if ((*it) == event.dut->getDUTid()) {
			continue;
		}
		transX = 0;
		transY = 0;
		//Translate if translations are found
		map<int, TransCalib>::iterator translations = m_translations.find(
				(*it));
		if (translations != m_translations.end()) {
			transX = (*translations).second.shiftX;
			transY = (*translations).second.shiftY;
		}
		//Loop over all tracks (Parameters?)
		for (int ii = 0; ii < t_nTrackParams; ii++) {
			// do only look at the TrackParams for the current DUT
			if (t_iden->at(ii) != (*it)) {
				continue;
			}
			bool foundMatch = false;	// for each trackParam reset foundMatch
			// the above line could be moved before the loop for clarity,
			// since this loop is terminated by a break; anyway in case of
			// foundMatch == true (see end of loop)
			const DUT* dut = config.getDut((*it));
			assert(dut != NULL);	// abort tbmon if not a valid DUT
			// skip if not first matched matchDUT
			// and other match was found (redundant? I say yes! because as soon as there is a match the trackNum will be changed and first will be false)
			// and that is not the same Track ID for current matchDUT
			if ((not first) and trackNum > -1
					and t_trackNum->at(ii) != trackNum) {
				continue;
			}
			// skip TrackParam if no matched matchDUT yet and chi2 too high
			if (first and t_chi2->at(ii) > chi2Min) {
				continue;
			}
			// multiply eutrack positions time 1000 to get um values
			double xx(t_posX->at(ii) * 1000.0 - transX);
			double yy(t_posY->at(ii) * 1000.0 - transY);
			//Loop over all hits, match to tracks.
			for (vector<PllHit*>::iterator hit = m_hitMap[(*it)].begin();
					hit != m_hitMap[(*it)].end(); hit++) {
				int col = (*hit)->col;
				int row = (*hit)->row;
				// skip if hit is outside RefLimitX in x direction
				if (fabs(col * dut->getPitchX() - xx)
						> dut->getRefLimitX() * dut->getPitchX()) {
					continue;
				}
				// skip if hit is outside RefLimitY in y direction
				if (fabs(row * dut->getPitchY() - yy)
						> dut->getRefLimitY() * dut->getPitchY()) {
					continue;
				}
				// if other matchDUT was already matched or
				// no other matchDUT yet and no match in this one yet: nMatch++
				// (this causes maximum nMatch == 1 for first matchDUT)
				// -> also trackNum will be the first found matching one
				if ((not first) or nMatch == 0)
					nMatch++;
				foundMatch = true;
				trackNum = t_trackNum->at(ii);	// save Track ID of match
				// if other matchDUT was already matched skip the rest
				// of the hits (this ensures max. increase of nMatch of 1
				// for each DUT)
				if (not first) {
					break;
				}
			}
			// if match found, skip all other trackParams (this skips further
			// tracks, in case there are any)
			if (foundMatch) {
				break;
			}
		}
		// if match was found, declare first = false
		// -> further matchDUTs are only scanned at the same Track ID
		if (nMatch > 0) {
			first = false;
		}
	}
	//trackNum = 0;
	//Extract trackParams
	for (int ii = 0; ii < t_nTrackParams; ii++) {
		// now only look at DUT this event belongs to
		if (t_iden->at(ii) != event.dut->getDUTid()) {
			continue;
		}
		// and only at the matched track
		if (t_trackNum->at(ii) != trackNum) {
			continue;
		}
		// This sets base and track flag to good if found
		// and takes the trackParameter at the position of the Event DUT
		event.setEvent(t_posX->at(ii) * 1000.0, t_posY->at(ii) * 1000.0,
				t_dxdz->at(ii), t_dydz->at(ii), // changed by BJD from 0, 0,
				t_chi2->at(ii), t_ndof->at(ii));
		
		// get rotation data for current DUT and save into event
		vector<int>::iterator it;
		it = std::find(gr_ID->begin(), gr_ID->end(), t_iden->at(ii));
		//cout << "Searching rotation data for plane " <<  t_iden->at(ii) << endl;
		//cout << "Position in rotation vector " << std::distance(gr_ID->begin(),it) << endl;
		if (getHasRotationTree() == true) {
			event.gr_RotXY = gr_RotXY->at(std::distance(gr_ID->begin(), it));
			event.gr_RotZX = gr_RotZX->at(std::distance(gr_ID->begin(), it));
			event.gr_RotZY = gr_RotZY->at(std::distance(gr_ID->begin(), it));
		};
		//cout << event.gr_RotXY << endl;
		break;
	}
	//Translate the track parameters
	map<int, TransCalib>::iterator translations = m_translations.find(
			event.dut->getDUTid());
	if (translations != m_translations.end()) {
		event.trackX -= (*translations).second.shiftX;
		event.trackY -= (*translations).second.shiftY;
		
	}
	// Set flags
	if (nMatch < m_nMatch) {
		event.fTrackMatchNeighbour = event::kBad;
		event.fTrack = event::kBad;
	} else {
		event.fTrackMatchNeighbour = event::kGood;
	}
	// if(event.ndof < 6){
	//   event.fTrack = event::kBad;
	// }
	//event.fTrackMatchNeighbour = event::kGood;
	//event.fTrack = event::kGood;
}

void EuBuildTrack::addTranslation(int iden, int firstRun, int lastRun,
		double shiftX, double shiftY) {
	TransCalib calib;
	calib.iden = iden;
	calib.firstRun = firstRun;
	calib.lastRun = lastRun;
	calib.shiftX = shiftX;
	calib.shiftY = shiftY;
	m_calibs.push_back(calib);
}

void EuBuildTrack::addTranslation(const char* filename, int iden, double addX, double addY) {
	ifstream file(filename);
	if (file.fail()) {
		cout
				<< "Failed to open file for the initialization of translator for iden: "
				<< iden << endl;
		return;
	}
	int runnumber;
	while (!file.eof()) {
		TransCalib calib;
		file >> runnumber >> calib.shiftX >> calib.shiftY;
		calib.shiftX += addX;
		calib.shiftY += addY;
		calib.iden = iden;
		calib.firstRun = runnumber;
		calib.lastRun = runnumber;
		m_calibs.push_back(calib);
	}
	file.close();
}

void EuBuildTrack::addMatchDUT(int iden) {
	m_idens.push_back(iden);
}
void EuBuildTrack::nMatches(int matches) {
	m_nMatch = matches;
}
