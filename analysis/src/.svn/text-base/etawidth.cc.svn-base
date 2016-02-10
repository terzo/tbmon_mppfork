#include "etawidth.h"

void EtaWidth::initHisto(vector<TH1D*> &histoc, vector<TH1D*>& histob,
		double halfWidth) {

	for (int ii = 0; ii < 110; ii++) {
		TH1D* c = new TH1D("", "", 1000, -halfWidth, halfWidth);
		TH1D* b = new TH1D("", "", 1000, -halfWidth, halfWidth);
		histoc.push_back(c);
		histob.push_back(b);
	}

}

void EtaWidth::init(TbConfig &config) {

	initHisto(residualsY2c, residualsY2b, 200);
	initHisto(residualsY3b, residualsY3c, 200);

	residualY22 = new TH1D("", "", 100, -50, 50);
	residualY = new TH1D("", "", 100, -50, 50);
	etaY = new TH1D("", "", 100, -50, 50);

}

double EtaWidth::borderCorrection(double ecorr, double width) {

	double shift = 0.5 * (1.0 - (width / 100.0));
	double scaled = (width / 100.0) * ecorr;
	return (scaled + shift);

}

double EtaWidth::centralCorrection(double ecorr, double width) {

	double shift(0.0);
	if (ecorr > 0.5) {
		shift = 1.0 - width / 100.0;
	}
	double scaled = (width / 100.0) * ecorr;
	return (scaled + shift);

}

void EtaWidth::event(const TbConfig &config, const Event &event) {

	// Standard track quality cuts
	if (event.fTrack != event::kGood) {
		return;
	}
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	// Require exactly one cluster in event, matched to track
	if (event.clusters.size() != 1) {
		return;
	}
	int matched = cluster::getMatched(event.clusters, event);
	if (matched == -1) {
		return;
	}
	// Cluster must have sufficient ToT
	const event::cluster_t cluster = event.clusters.at(matched);
	if (cluster::getSumTot(cluster) < 3) {
		return;
	}

	double ecorr = cluster::getChargeWeightedRow(cluster);
	bool accept = true;
	if (fabs(event.trackY - event.dut->getPitchY() * ecorr) > 50.0)
		accept = false;
	double floored = floor(ecorr);
	ecorr = ecorr - floored;
	if (accept)
		etaY->Fill(
				event.trackY
						- cluster::getEtaCorrectedY(cluster, event, event.dut));
	if (cluster::spannedRows(cluster) == 2) {
		if (accept) {
			residualY->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + borderCorrection(ecorr, 38)));
			residualY22->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + borderCorrection(ecorr, 99)));
		}
		for (int ii = 0; ii < 110; ii++) {
			residualsY2c.at(ii)->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + centralCorrection(ecorr, ii)));
			residualsY2b.at(ii)->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + borderCorrection(ecorr, ii)));
		}
	} else if (cluster::spannedRows(cluster) == 3) {
		if (accept) {
			residualY->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + centralCorrection(ecorr, 80)));
			residualY22->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + centralCorrection(ecorr, 48)));
		}
		for (int ii = 0; ii < 110; ii++) {
			residualsY3c.at(ii)->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + centralCorrection(ecorr, ii)));
			residualsY3b.at(ii)->Fill(
					event.trackY
							- event.dut->getPitchY()
									* (floored + borderCorrection(ecorr, ii)));
		}
	} else {
		if (accept) {
			residualY22->Fill(
					event.trackY - cluster::getUnWeightedY(cluster, event));
			residualY->Fill(
					event.trackY - cluster::getUnWeightedY(cluster, event));
		}
	}

}

void EtaWidth::getStats(const TbConfig& config, vector<TH1D*> histos,
		const char* name) {

	TH1D* mean = new TH1D("", "", 110, 0, 110);
	mean->GetXaxis()->SetTitle("band width");
	mean->GetYaxis()->SetTitle("mean");
	TH1D* rms = new TH1D("", "", 110, 0, 110);
	rms->GetXaxis()->SetTitle("band width");
	rms->GetYaxis()->SetTitle("rms");
	for (int ii = 0; ii < 110; ii++) {
		mean->SetBinContent(ii + 1, histos.at(ii)->GetMean());
		rms->SetBinContent(ii + 1, histos.at(ii)->GetRMS());
	}
	char fname[400];
	sprintf(fname, "%s-rms", name);
	config.setTitle(this->name, fname, rms);
	config.drawToFile(this->name, fname, rms);
	sprintf(fname, "%s-mean", name);
	config.setTitle(this->name, fname, mean);
	config.drawToFile(this->name, fname, mean);

}

void EtaWidth::finalize(const TbConfig &config) {

	getStats(config, residualsY2c, "center2");
	getStats(config, residualsY2b, "border2");
	getStats(config, residualsY3c, "center3");
	getStats(config, residualsY3b, "border3");
	config.setTitle(this->name, "residual", residualY);
	config.drawToFile(this->name, "residual", residualY);
	TBALOG(kINFO) << "New RMS(0):" << residualY->GetRMS() << ","
			<< residualY->GetEntries() << endl;
	TBALOG(kINFO) << "New RMS(22):" << residualY22->GetRMS() << ","
			<< residualY22->GetEntries() << endl;
	TBALOG(kINFO) << "Eta RMS:" << etaY->GetRMS() << "," << etaY->GetEntries()
			<< endl;

}
