#include "getetacorr.h"

void GetEtaCorr::init(TbConfig &config) {

	h_chargeWeightedX = new TH1D("", ";Charge weighted cluster x position;Clusters", 101, 0.0, 1.01);
	h_chargeWeightedY = new TH1D("", ";Charge weighted cluster y position;Clusters", 101, 0.0, 1.01);
	h_etacorrsX = new TH1D("", ";Eta correction in x;Entries", 101, 0.0, 1.01);
	h_etacorrsY = new TH1D("", ";Eta correction in y;Entries", 101, 0.0, 1.01);

}

void GetEtaCorr::event(const TbConfig &config, const Event &event) {

	// Only one cluster, size of 1 or 2, non-zero ToT,
	// matched to a track passing through a good sensor region
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	if (event.fTrack != event::kGood) {
		return;
	}
	if (event.clusters.size() != 1) {
		return;
	}
	int iMatch = cluster::getMatched(event.clusters, event);
	if (iMatch != 0) {
		return;
	}
	if (event.clusters.at(iMatch).size() > 2) {
		return;
	}
	int sumTot = cluster::getSumTot(event.clusters.at(iMatch));
	if (sumTot == 0) {
		return;
	}

	double qmeanX = cluster::getChargeWeightedCol(event.clusters.at(iMatch));
	double qmeanY = cluster::getChargeWeightedRow(event.clusters.at(iMatch));
	h_chargeWeightedX->Fill(qmeanX - floor(qmeanX));
	h_chargeWeightedY->Fill(qmeanY - floor(qmeanY));

}

void GetEtaCorr::finalize(const TbConfig &config) {

	char* streamName = config.getOutStreamName(name, "etacalib");
	ofstream stream;
	stream.open(streamName);
	int lastrun = (config.lastRun == -1) ? config.firstRun : config.lastRun;
	stream << config.firstRun << " " << lastrun << " " << endl;
	printEtaCalib(h_chargeWeightedY, h_etacorrsY, stream);
	printEtaCalib(h_chargeWeightedX, h_etacorrsX, stream);

	config.drawAndSave(name, "chargeWeigtedY", h_chargeWeightedY);
	config.drawAndSave(name, "etacorrsY", h_etacorrsY);
	config.drawAndSave(name, "chargeWeigtedX", h_chargeWeightedX);
	config.drawAndSave(name, "etacorrsX", h_etacorrsX);

}

void GetEtaCorr::printEtaCalib(TH1D* histo, TH1D* etacorrs, ofstream &stream) {

	double singleHits = histo->GetBinContent(1);
	histo->SetBinContent(1, singleHits * 0.5);
	histo->SetBinContent(histo->GetNbinsX(), singleHits * 0.5);
	histo->Scale(1.0 / histo->Integral());
	double integral = 0; //histo->GetBinContent(1);
	for (int bin = 1; bin < histo->GetNbinsX(); bin++) {
		integral += histo->GetBinContent(bin);
		etacorrs->SetBinContent(bin, integral);
		stream << integral << " ";
	}
	//Remove single hits for plotting
	histo->SetBinContent(1, 0.0);
	histo->SetBinContent(histo->GetNbinsX(), 0.0);

}

