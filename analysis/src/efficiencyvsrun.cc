#include "efficiencyvsrun.h"

void EfficiencyVsRun::init(TbConfig & config) {

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	//
	// initialize maps with n(tracks), n(any hits), n(matched hits) vs run
	//
	//int nRuns = config.runList.size();
	for (std::vector<int>::const_iterator currentRun = config.runList.begin();
			currentRun != config.runList.end(); currentRun++) {

		ntracks[*currentRun] = 0;
		nhits[*currentRun] = 0;
		nmatched[*currentRun] = 0;
	}

	//Get -P:key val cmd line arguments
	char key[100];
	sprintf(key, "efficiencyvsrun_effPlotMinY_%u", this->iden);
	effPlotMinY = config.cmdLineExtras_argGetter(key, 0.9,
			"Minimum plotted value in efficiency plots");
}

void EfficiencyVsRun::event(const TbConfig &config, const Event &event) {

	TBALOG(kDEBUG3) << "Started analysis class event() command" << endl;
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		TBALOG (kDEBUG3) << "Event not passed all cuts" << endl;
		return;
	}
	// Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		TBALOG (kDEBUG3) << "Clusters not present in event" << endl;
		return;
	}
	// Are we in the central region on the chip?
	if (event.fTrackCentralRegion != event::kGood) {
		return;
	}
	// Are we in a good region of the chip?
	if (event.fTrackRegion != event::kGood) {
		return;
	}

	ntracks[config.currentRun]++;

	TBALOG(kDEBUG2) << "Event: " << config.currentEntry << " Run: "
			<< config.currentRun << endl;
	TBALOG(kDEBUG2) << "Accepted track" << endl;

	// Are there any hits in the DUT
	if (event.hits.size() > 0) {
		nhits[config.currentRun]++;
		TBALOG(kDEBUG2) << "Found hit in DUT!" << endl;
	}

	// Look for matched clusters
	if (cluster::getMatched(event.clusters, event) != -1) {
		nmatched[config.currentRun]++;
		TBALOG(kDEBUG2) << "Found matching cluster!" << endl;
	}
}

void EfficiencyVsRun::finalize(const TbConfig &config) {

	//
	// initialize graphs (MOVE TO END)
	//
	int nRuns = ntracks.size();

	double* runs = new double[nRuns];
	double* effall = new double[nRuns];
	double* e_effall = new double[nRuns];
	double* effmatched = new double[nRuns];
	double* e_effmatched = new double[nRuns];

	histo_EffAllVsRun = new TH1D("", ";Run Number;(All) Efficiency %",
			config.runList.at(nRuns - 1) - config.runList.at(0),
			config.runList.at(0), config.runList.at(nRuns - 1) + 1);
	histo_EffMatchedVsRun = new TH1D("", ";Run Number;(Matched) Efficiency %",
			config.runList.at(nRuns - 1) - config.runList.at(0),
			config.runList.at(0), config.runList.at(nRuns - 1) + 1);

	for (int irun = 0; irun < nRuns; irun++) {

		int run = config.runList.at(irun);
		runs[irun] = (double) run;
		double trk = (double) (ntracks[run]);
		double hit = (double) (nhits[run]);
		double matched = (double) (nmatched[run]);
		if (trk > 0) {
			effall[irun] = hit / trk;
			if (effall[irun] < 1.0)
				e_effall[irun] = sqrt(effall[irun] * (1. - effall[irun]) / trk);
			else
				e_effall[irun] = 0.0;

			effmatched[irun] = matched / trk;
			if (effmatched[irun] < 1.0)
				e_effmatched[irun] = sqrt(
						effmatched[irun] * (1. - effmatched[irun]) / trk);
			else
				e_effmatched[irun] = 0.0;

		} else {
			effall[irun] = 0.0;
			e_effall[irun] = 0.0;

			effmatched[irun] = 0.0;
			e_effmatched[irun] = 0.0;
		}
		histo_EffAllVsRun->SetBinContent(run - config.runList.at(0) + 1,
				effall[irun]);
		histo_EffAllVsRun->SetBinError(run - config.runList.at(0) + 1,
				e_effall[irun]);
		histo_EffMatchedVsRun->SetBinContent(run - config.runList.at(0) + 1,
				effmatched[irun]);
		histo_EffMatchedVsRun->SetBinError(run - config.runList.at(0) + 1,
				e_effmatched[irun]);
	}

	config.drawAndSave(this->name, (const char*) "effAllVsRun",
			histo_EffAllVsRun);
	config.drawAndSave(this->name, (const char*) "effMatchedVsRun",
			histo_EffMatchedVsRun);

	h_EffAllVsRun = new TGraphErrors(nRuns, runs, effall, 0, e_effall);
	h_EffMatchedVsRun = new TGraphErrors(nRuns, runs, effmatched, 0,
			e_effmatched);
	drawAndSave(config, name, "effAllVsRun", h_EffAllVsRun);
	drawAndSave(config, name, "effMatchedVsRun", h_EffMatchedVsRun);

	delete[] runs;
	delete[] effall;
	delete[] e_effall;
	delete[] effmatched;
	delete[] e_effmatched;

	/*
	 TBALOG(kINFO) << "Tracks:             " << ntracks << endl;
	 TBALOG(kINFO) << "Hits:               " << anyhit << endl;
	 TBALOG(kINFO) << "Tracks w/ hit:      " << nhits << endl;
	 TBALOG(kINFO) << "EfficiencyVsRun (match): " << double(nhits)/double(ntracks) << endl;
	 TBALOG(kINFO) << "EfficiencyVsRun (any):   " << double(anyhit)/double(ntracks) << endl;
	 */
}

void EfficiencyVsRun::drawAndSave(const TbConfig& config,
		const char* analysisName, const char* histoName, TGraphErrors* gr) {

	char* fileName = config.buildHistName(analysisName, histoName);
	const DUT* dut = config.getDut(this->iden);

	TCanvas* c_eff = new TCanvas("c_eff", "c_eff");
	c_eff->cd();

	gr->Draw("AP");
	gr->GetXaxis()->SetTitle("Run number");
	gr->GetYaxis()->SetTitle("Efficiency");
	//gr->GetYaxis()->SetRangeUser(0.90,1.01);
	gr->GetYaxis()->SetRangeUser(effPlotMinY, 1.01);
	c_eff->SaveAs(fileName);

	delete c_eff;
}
