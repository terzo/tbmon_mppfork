#include "checkalign.h"

void CheckAlign::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	h_hitxVresx = new TProfile("", ";Hit Position [col];X Residual [#mum]",
			nCols, 0, nCols);
	h_hityVresx = new TProfile("", ";Hit Position [row];X Residual [#mum]",
			nRows, 0, nRows);
	h_hitxVresy = new TProfile("", ";Hit Position [col];Y Residual [#mum]",
			nCols, 0, nCols);
	h_hityVresy = new TProfile("", ";Hit Position [row];Y Residual [#mum]",
			nRows, 0, nRows);

	doCuts = true;

}

void CheckAlign::initRun(const TbConfig &config) {

	h_resX = new TH1D("", ";X Residual [#mum]", 2000, -1000, 1000);
	h_resY = new TH1D("", ";Y Residual [#mum]", 2000, -1000, 1000);

}

void CheckAlign::event(const TbConfig &config, const Event& event) {

	if (doCuts == true) {
		// Standard track quality cuts
		if (event.fTrack != event::kGood) {
			return;
		}
		if (event.fTrackRegion != event::kGood) {
			return;
		}
	}

	for (vector<event::cluster_t>::const_iterator iClu = event.clusters.begin();
			iClu != event.clusters.end(); iClu++) {
		if (event.clusters.size() == 1 and cluster::getSumTot(*iClu) >= 3) {
			h_resX->Fill(
					event.trackX
							- cluster::getEtaCorrectedX(*iClu, event,
									event.dut));
			h_resY->Fill(
					event.trackY
							- cluster::getEtaCorrectedY(*iClu, event,
									event.dut));
		}
		if (doCuts == true) {
			if (cluster::isMatch(*iClu, event) == 0) {
				return;
			}
		}
		if (cluster::getSumTot(*iClu) == 0) {
			return;
		}

		double hitX = cluster::getChargeWeightedCol(*iClu);
		double hitY = cluster::getChargeWeightedRow(*iClu);
		if (doCuts == true) {
			if (hitX < event.dut->getSkipCols()
					or hitX
							> (event.dut->getNcols() - event.dut->getSkipCols())) {
				return;
			}
			if (hitY < event.dut->getSkipRows()
					or hitY
							> (event.dut->getNrows() - event.dut->getSkipRows())) {
				return;
			}
		}

		double resX = event.trackX - (hitX * event.dut->getPitchX());
		double resY = event.trackY - (hitY * event.dut->getPitchY());
		if (doCuts == false or (doCuts == true and event.ndof >= 2)) {
			h_hitxVresx->Fill(hitX, resX);
			h_hityVresx->Fill(hitY, resX);
			h_hitxVresy->Fill(hitX, resY);
			h_hityVresy->Fill(hitY, resY);
		}
	}
	/* BJD
	 // Below is the old (cleaned-up) version! Cuts required, only 1 matching cluster used
	 // Above gives the option of turning off cuts, and also checking multiple clusters,
	 // which effectively combines the checkalign and checkalignsimp classes
	 if(event.fTrack != event::kGood) { return; }
	 if( event.fTrackRegion != event::kGood) { return; }

	 int matched = cluster::getMatched(event.clusters, event);
	 if( matched == -1 ) { return; }
	 if( cluster::getSumTot( event.clusters.at(matched)) == 0) { return; }

	 double hitX = cluster::getChargeWeightedCol( event.clusters.at(matched));
	 double hitY = cluster::getChargeWeightedRow( event.clusters.at(matched));
	 //if( hitX < 2 or hitX > 15) { return; }
	 //if( hitY < 16 or hitY > 144) { return; } // BJD
	 if ( hitX < event.dut->getSkipCols() or hitX > (event.dut->getNcols() - event.dut->getSkipCols()) ) { return; }
	 if ( hitY < event.dut->getSkipRows() or hitY > (event.dut->getNrows() - event.dut->getSkipRows()) ) { return; }

	 double resX = event.trackX - (hitX * event.dut->getPitchX());
	 double resY = event.trackY - (hitY * event.dut->getPitchY());
	 if(event.ndof >= 2){
	 h_hitxVresx->Fill(hitX, resX);
	 h_hityVresx->Fill(hitY, resX);
	 h_hitxVresy->Fill(hitX, resY);
	 h_hityVresy->Fill(hitY, resY);
	 }
	 */
}

void CheckAlign::finalizeRun(const TbConfig &config) {

	if (fabs(h_resX->GetMean()) > 10 or fabs(h_resY->GetMean()) > 10) {
		badRuns.push_back(config.currentRun);
		char title[500];
		sprintf(title, "run%iX", config.currentRun);
		config.saveToFile(this->name, title, h_resX);
		sprintf(title, "run%iY", config.currentRun);
		config.saveToFile(this->name, title, h_resY);
	}

	deltaX.push_back(h_resX->GetMean());
	deltaY.push_back(h_resY->GetMean());
	run.push_back(config.currentRun);

	delete h_resY, h_resX;

}

void CheckAlign::finalize(const TbConfig &config) {

	config.drawAndSave(this->name, "hitXresX", h_hitxVresx);
	config.drawAndSave(this->name, "hitYresX", h_hityVresx);
	config.drawAndSave(this->name, "hitXresY", h_hitxVresy);
	config.drawAndSave(this->name, "hitYresY", h_hityVresy);

	//Fit a pol 1 to get a measure for the misalignment and the tilt of the sample
	//Prepare the functions and fit
	TF1 *pol1_xx = new TF1("pol1_xx", "pol1", 3, 15);
	h_hitxVresx->Fit("pol1_xx", "+ QR");
	TF1 *pol1_yx = new TF1("pol1_yx", "pol1", 15, 135);
	h_hityVresx->Fit("pol1_yx", "+ QR");
	TF1 *pol1_xy = new TF1("pol1_xy", "pol1", 3, 15);
	h_hitxVresy->Fit("pol1_xy", "+ QR");
	TF1 *pol1_yy = new TF1("pol1_yy", "pol1", 15, 135);
	h_hityVresy->Fit("pol1_yy", "+ QR");
	//Pipe out the plots
	config.drawAndSave(this->name, "hitXresX_fitted", h_hitxVresx);
	config.drawAndSave(this->name, "hitYresX_fitted", h_hityVresx);
	config.drawAndSave(this->name, "hitXresY_fitted", h_hitxVresy);
	config.drawAndSave(this->name, "hitYresY_fitted", h_hityVresy);
	TBALOG(kINFO) << "hitXresX Offset: " << pol1_xx->GetParameter(0) << " +- "
			<< pol1_xx->GetParError(0) << " µm, tilt: "
			<< pol1_xx->GetParameter(1) << " +- " << pol1_xx->GetParError(1)
			<< "\n";
	TBALOG(kINFO) << "hitYresX Offset: " << pol1_yx->GetParameter(0) << " +- "
			<< pol1_yx->GetParError(0) << " µm, tilt: "
			<< pol1_yx->GetParameter(1) << " +- " << pol1_yx->GetParError(1)
			<< "\n";
	TBALOG(kINFO) << "hitXresY Offset: " << pol1_xy->GetParameter(0) << " +- "
			<< pol1_xy->GetParError(0) << " µm, tilt: "
			<< pol1_xy->GetParameter(1) << " +- " << pol1_xy->GetParError(1)
			<< "\n";
	TBALOG(kINFO) << "hitYresY Offset: " << pol1_yy->GetParameter(0) << " +- "
			<< pol1_yy->GetParError(0) << " µm, tilt: "
			<< pol1_yy->GetParameter(1) << " +- " << pol1_yy->GetParError(1)
			<< "\n";

	if (badRuns.size() == 0) {
		TBALOG(kINFO) << "Found 0 suspicious runs!\n";
	} else if (badRuns.size() > 0) {
		TBALOG(kINFO) << "Found " << badRuns.size() << " suspicious runs:\n";
		for (int ii = 0; ii < badRuns.size(); ii++) {
			cout << badRuns.at(ii) << endl;
			cout << "delta X [um] = " << deltaX.at(ii) << endl;
			cout << "delta Y [um] = " << deltaY.at(ii) << endl;
		}
	}

	char* streamName = config.getOutStreamName(name, "translations");
	ofstream stream;
	stream.open(streamName);

	for (int ii = 0; ii < run.size(); ii++) {
		stream << run.at(ii) << " " << deltaX.at(ii) << " " << deltaY.at(ii)
				<< endl;
	};

	stream.close();

}
