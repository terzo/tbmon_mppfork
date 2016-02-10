// Log 
// 2011_Aug_05: Modified by Sh.Tsiskaridze
// 2011_Sep_19: Modified by Sh.Tsiskaridze; Right edge of the sensors are added; 

#include "edgeefficiency.h"

void EdgeEfficiency::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();
	double epitchX = dut->getePitchX();
	double epitchY = dut->getePitchY();
	int skipCols = dut->getSkipCols();
	int skipRows = dut->getSkipRows();
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	nTracks = 0;
	nTracksEdge = 0;
	nHitsEdge = 0;

	h_EEAnyHitSensorMap = new TH2D("", "Any Edge Hit Map;Column;Row", nCols, 0,
			nCols, nRows, 0, nRows);
	h_EEHitSensorMap = new TH2D("", "Hit Map Edge Tracks;Column;Row", nCols + 4,
			-2, nCols + 2, nRows, 0, nRows);
	h_EETrackSensorMap = new TH2D("", "Track Map Edge Tracks;Column;Row",
			nCols + 4, -2, nCols + 2, nRows, 0, nRows);

	h_EESensorMap = new TH2D("", "Hit Map Edge Tracks;Column;Row", nCols + 4,
			-2, nCols + 2, nRows, 0, nRows);

	h_EEHitPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EEHitLeftPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EEHitRightPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);

	h_EETrackPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EETrackLeftPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EETrackRightPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
}

void EdgeEfficiency::event(const TbConfig &config, const Event &event) {

	double pitchX = event.dut->getPitchX();
	double pitchY = event.dut->getPitchY();
	double epitchX = event.dut->getePitchX();
	double epitchY = event.dut->getePitchY();
	int skipCols = event.dut->getSkipCols();
	int skipRows = event.dut->getSkipRows();
	int nCols = event.dut->getNcols();
	int nRows = event.dut->getNrows();

	//Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		TBALOG (kDEBUG3) << "Event did not pass all cuts" << endl;
		return;
	}
	//Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		TBALOG (kDEBUG3) << "Clusters not present in event" << endl;
		return;
	}

	nTracks++;

	TBALOG(kDEBUG2) << "Event: " << config.currentEntry << " Run: "
			<< config.currentRun << endl;
	TBALOG(kDEBUG2) << "Good track." << endl;

	// Are we *not* in the central region on the chip? 
	if (event.fTrackCentralRegion == event::kGood) {
		TBALOG (kDEBUG3) << "Track in central region of chip" << endl;
		return;
	}

	TBALOG(kDEBUG2) << "Edge track, col = "
			<< tbutils::getCol(event.trackX, event) << " row = "
			<< tbutils::getRow(event.trackY, event) << endl;

	// The edge in Y (ganged- and inter-ganged pixels) should still be ignored
	if ((event.trackY + 0.5 * pitchY) < (pitchY * skipRows)) {
		return;
	}
	if ((event.trackY + 0.5 * pitchY) > (pitchY * nRows - skipRows)) {
		return;
	}

	// Make sure we are not too far off the sensor in the long (X) direction for edge studies
	// (for the X edge we want to include half a pixel out)  
	// select tracks for LEFT edge of the sensor
	if ((-1.0 * epitchX < event.trackX) && (event.trackX < pitchX)) {
		h_EETrackPixelMap->Fill(event.trackX + epitchX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));
		h_EETrackLeftPixelMap->Fill(event.trackX + epitchX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));

		// Look for matched clusters only
		if (cluster::getMatched(event.clusters, event) != -1) {
			TBALOG(kDEBUG2) << "Found matching cluster!" << endl;
			h_EEHitPixelMap->Fill(event.trackX + epitchX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
			h_EEHitLeftPixelMap->Fill(event.trackX + epitchX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
		}
	} else
	// select tracks for RIGHT edge of the sensor
	if (((nCols - 2) * pitchX < event.trackX)
			&& (event.trackX < (nCols - 1) * pitchX + epitchX)) {
		h_EETrackPixelMap->Fill(((nCols - 1) * pitchX + epitchX) - event.trackX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));
		h_EETrackRightPixelMap->Fill(event.trackX - (nCols - 2) * pitchX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));

		// Look for matched clusters only
		if (cluster::getMatched(event.clusters, event) != -1) {
			TBALOG(kDEBUG2) << "Found matching cluster!" << endl;
			h_EEHitPixelMap->Fill(
					((nCols - 1) * pitchX + epitchX) - event.trackX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
			h_EEHitRightPixelMap->Fill(event.trackX - (nCols - 2) * pitchX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
		}
	} else {
		return;
	}

	// Fill Sensor Maps using both sides of the sensor
	// Count the tracks
	nTracksEdge++;

	TBALOG(kDEBUG2) << "Event: " << config.currentEntry << " Run: "
			<< config.currentRun << endl;
	TBALOG(kDEBUG2) << "Accepted track." << endl;

	// plot the position of *any* single hits
	for (vector<PllHit*>::const_iterator i = event.hits.begin();
			i != event.hits.end(); ++i) {
		if ((*i)->iden != this->iden)
			continue;
		h_EEAnyHitSensorMap->Fill((*i)->col, (*i)->row);
	}

	float col = (event.trackX + 0.5 * pitchX) / pitchX;
	float row = (event.trackY + 0.5 * pitchY) / pitchY;

	// Fill the Edge Efficiency Sensor Map
	h_EETrackSensorMap->Fill(col, row);

	// Are there any hits in the DUT?
	if (event.hits.size() > 0) {
		TBALOG(kDEBUG2) << "Found hit in DUT!" << endl;
	}

	// Look for matched clusters
	if (cluster::getMatched(event.clusters, event) == -1) {
		return;
	}
	TBALOG(kDEBUG2) << "Found matching cluster!" << endl;

	// count the hits 
	nHitsEdge++;
	h_EEHitSensorMap->Fill(col, row);
	h_EESensorMap->SetStats(false);
	h_EESensorMap->Divide(h_EEHitSensorMap, h_EETrackSensorMap, 1.0, 1.0, "B");

}

void EdgeEfficiency::finalize(const TbConfig &config) {
	char* fileName;
	TBALOG(kINFO) << " Tracks:             " << nTracks << endl;
	TBALOG(kINFO) << " Edge tracks:        " << nTracksEdge << endl;
	TBALOG(kINFO) << " Edge tracks w/ hit: " << nHitsEdge << endl;

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double epitchX = dut->getePitchX();

	//create Sensor Map files

	fileName = config.buildHistName(name, "edgeAnyHitSensorMap");
	ApplySensorMapStyle(config, h_EEAnyHitSensorMap, "Edge Any Hit Sensor Map",
			fileName, "edgeAnyHitSensorMap");

	fileName = config.buildHistName(name, "edgeHitSensorMap");
	ApplySensorMapStyle(config, h_EEHitSensorMap, "Edge Hit Sensor Map",
			fileName, "edgeHitSensorMap");

	fileName = config.buildHistName(name, "edgeTrackSensorMap");
	ApplySensorMapStyle(config, h_EETrackSensorMap, "Edge Track Sensor Map",
			fileName, "edgeTrackSensorMap");

	fileName = config.buildHistName(name, "edgeEffSensorMap");
	ApplySensorMapStyle(config, h_EESensorMap, "Edge Efficiency Sensor Map",
			fileName, "edgeEffSensorMap");

	//create Pixel Map files
	fileName = config.buildHistName(name, "edgeTrackPixelMap");
	ApplyPixelMapStyle(config, h_EETrackPixelMap, "Edge Track Pixel Map",
			fileName, "edgeTrackPixelMap");

	fileName = config.buildHistName(name, "edgeTrackLeftPixelMap");
	ApplyPixelMapStyle(config, h_EETrackLeftPixelMap,
			"Left Edge Tracks Pixel Map", fileName, "edgeTrackLeftPixelMap");

	fileName = config.buildHistName(name, "edgeTrackPixelMap");
	ApplyPixelMapStyle(config, h_EETrackRightPixelMap,
			"Right Edge Tracks Pixel Map", fileName, "edgeTrackRightPixelMap");

	fileName = config.buildHistName(name, "edgeHitPixelMap");
	ApplyPixelMapStyle(config, h_EEHitPixelMap, "Edge Hit Pixel Map", fileName,
			"edgeHitPixelMap");

	fileName = config.buildHistName(name, "edgeHitLeftPixelMap");
	ApplyPixelMapStyle(config, h_EEHitLeftPixelMap, "Left Edge Hit Pixel Map",
			fileName, "edgeHitLeftPixelMap");

	fileName = config.buildHistName(name, "edgeHitPixelMap");
	ApplyPixelMapStyle(config, h_EEHitRightPixelMap, "Right Edge Hit Pixel Map",
			fileName, "edgeHitRightPixelMap");

	// build Edge Efficiency 

	h_EETrackPixelMap->Rebin2D(5, 2);
	h_EEHitPixelMap->Rebin2D(5, 2);

	h_EETrackLeftPixelMap->Rebin2D(5, 2);
	h_EEHitLeftPixelMap->Rebin2D(5, 2);

	h_EETrackRightPixelMap->Rebin2D(5, 2);
	h_EEHitRightPixelMap->Rebin2D(5, 2);

	TH2D* h_EEPixelMap = (TH2D*) h_EEHitPixelMap->Clone(
			TString::Format("edgeEffPixelMap"));
	h_EEPixelMap->Divide(h_EEHitPixelMap, h_EETrackPixelMap, 1., 1., "B");
	fileName = config.buildHistName(name, "edgeEffPixelMap");
	ApplyPixelMapStyle(config, h_EEPixelMap, "Edge Efficiency Pixel Map",
			fileName, "edgeEffPixelMap");

	TH2D* h_EELeftPixelMap = (TH2D*) h_EEHitLeftPixelMap->Clone(
			TString::Format("edgeEffLeftPixelMap"));
	h_EELeftPixelMap->Divide(h_EEHitLeftPixelMap, h_EETrackLeftPixelMap, 1., 1.,
			"B");
	fileName = config.buildHistName(name, "edgeEffLeftPixelMap");
	ApplyPixelMapStyle(config, h_EELeftPixelMap,
			"Left Edge Efficiency Pixel Map", fileName, "edgeEffLeftPixelMap");

	TH2D* h_EERightPixelMap = (TH2D*) h_EEHitRightPixelMap->Clone(
			TString::Format("edgeEffRightPixelMap"));
	h_EERightPixelMap->Divide(h_EEHitRightPixelMap, h_EETrackRightPixelMap, 1.,
			1., "B");
	fileName = config.buildHistName(name, "edgeEffRightPixelMap");
	ApplyPixelMapStyle(config, h_EERightPixelMap,
			"Right Edge Efficiency Pixel Map", fileName,
			"edgeEffRightPixelMap");

	// both sides
	char projname[500];
	sprintf(projname, "trackEdgePixelMapX-%d", this->iden);
	TH1D* h_EETrackPixelMapX = (TH1D*) h_EETrackPixelMap->ProjectionX(projname,
			-1, 999999, "E");

	sprintf(projname, "hitEdgePixelMapX-%d", this->iden);
	TH1D* h_EEHitPixelMapX = (TH1D*) h_EEHitPixelMap->ProjectionX(projname, -1,
			999999, "E");
	TGraphAsymmErrors* tg_EEPixelMapX = new TGraphAsymmErrors(h_EEHitPixelMapX,
			h_EETrackPixelMapX, "cl=0.683 b(1,1) mode");

	drawToFileEdgeEff(config, name, "edgeEff", tg_EEPixelMapX, 10,
			(pitchX + epitchX) / 2);

	//Left side
	char projnameL[500];
	sprintf(projnameL, "trackEdgeLeftPixelMapX-%d", this->iden);
	TH1D* h_EETrackLeftPixelMapX = (TH1D*) h_EETrackLeftPixelMap->ProjectionX(
			projnameL, -1, 999999, "E");

	sprintf(projnameL, "hitEdgeLeftPixelMapX-%d", this->iden);
	TH1D* h_EEHitLeftPixelMapX = (TH1D*) h_EEHitLeftPixelMap->ProjectionX(
			projnameL, -1, 999999, "E");
	TGraphAsymmErrors* tg_EELeftPixelMapX = new TGraphAsymmErrors(
			h_EEHitLeftPixelMapX, h_EETrackLeftPixelMapX,
			"cl=0.683 b(1,1) mode");

	drawToFileEdgeEff(config, name, "edgeEffLeft", tg_EELeftPixelMapX, 10,
			(pitchX + epitchX) / 2);

	//Right side
	char projnameR[500];
	sprintf(projnameR, "trackEdgeRightPixelMapX-%d", this->iden);
	TH1D* h_EETrackRightPixelMapX = (TH1D*) h_EETrackRightPixelMap->ProjectionX(
			projnameR, -1, 999999, "E");

	sprintf(projname, "hitEdgeRightPixelMapX-%d", this->iden);
	TH1D* h_EEHitRightPixelMapX = (TH1D*) h_EEHitRightPixelMap->ProjectionX(
			projname, -1, 999999, "E");
	TGraphAsymmErrors* tg_EERightPixelMapX = new TGraphAsymmErrors(
			h_EEHitRightPixelMapX, h_EETrackRightPixelMapX,
			"cl=0.683 b(1,1) mode");

	drawToFileEdgeEff(config, name, "edgeEffRight", tg_EERightPixelMapX,
			(pitchX + epitchX) / 2, pitchX + epitchX - 10);

	// delete histos
	delete fileName;
	delete h_EETrackSensorMap;
	delete h_EEAnyHitSensorMap;
	delete h_EEHitSensorMap;
	delete h_EESensorMap;

	delete h_EETrackPixelMap;
	delete h_EETrackLeftPixelMap;
	delete h_EETrackRightPixelMap;

	delete h_EEHitPixelMap;
	delete h_EEHitLeftPixelMap;
	delete h_EEHitRightPixelMap;

	delete h_EEPixelMap;
	delete h_EELeftPixelMap;
	delete h_EERightPixelMap;

	delete h_EETrackPixelMapX;
	delete h_EETrackLeftPixelMapX;
	delete h_EETrackRightPixelMapX;

	delete tg_EEPixelMapX;
	delete tg_EELeftPixelMapX;
	delete tg_EERightPixelMapX;

}

void EdgeEfficiency::drawToFileEdgeEff(const TbConfig& config,
		const char* analysisName, const char* histoName, TObject* histo,
		double fitLeft, double fitRight) const {

	const DUT* dut = config.getDut(this->iden);
	char* fileName = config.buildHistName(analysisName, histoName);
	double pitchX = dut->getPitchX();
	double epitchX = dut->getePitchX();

	// Make sure both a graph and hist can be used in the fit
	TGraphAsymmErrors* gr_histo = 0;
	TH1D* h_histo = 0;
	bool isGraph = false;
	if (histo->InheritsFrom("TGraph")) {
		isGraph = true;
		gr_histo = (TGraphAsymmErrors*) histo;
	} else if (histo->InheritsFrom("TH1D")) {
		h_histo = (TH1D*) histo;
	} else {
		TBALOG(kERROR)
				<< "You must have either a TH1D or TGraph for the edge fit!"
				<< endl;
		exit(1);
	}
	assert(tbutils::XOR(gr_histo, h_histo));

	TCanvas* canvas = new TCanvas("", "", 1200, 300);
	canvas->cd();
	gStyle->SetPadBottomMargin(0.15);
	gROOT->ForceStyle();
	TString title = ";Long pixel [#mum ]; Edge efficiency";
	TH1* h_plot = 0;
	if (isGraph) {
		gr_histo->SetLineWidth((Width_t) 1.3); // BJD: Really?!
		gr_histo->SetMarkerStyle(20);
		gr_histo->SetMarkerSize(0.7);
		gr_histo->Draw("AP");
		gr_histo->SetTitle(title);
		gr_histo->GetXaxis()->SetLabelSize(0.07);
		gr_histo->GetYaxis()->SetLabelSize(0.07);
		gr_histo->GetYaxis()->SetTitleSize(0.08);
		gr_histo->GetXaxis()->SetTitleSize(0.08);
		gr_histo->GetYaxis()->SetTitleOffset(0.4);
		gr_histo->GetXaxis()->SetTitleOffset(0.9);
		gr_histo->GetXaxis()->SetRangeUser(0., pitchX + epitchX);
		h_plot = gr_histo->GetHistogram();
	} else {
		h_histo->SetLineWidth((Width_t) 1.3); // BJD: Really?!
		h_histo->SetTitle(title);
		h_histo->GetYaxis()->SetTitleSize(0.09);
		h_histo->GetYaxis()->SetTitleOffset(0.04);
		h_histo->Draw();
		h_plot = h_histo;
	}

	// Draw an area covering the next to last pixel for visualization
	TH1D* h_innerPixel = (TH1D*) h_plot->Clone("innerPixel");
	for (int iBin = 1; iBin != h_innerPixel->GetNbinsX() + 1; ++iBin) {
		if ((iBin <= h_plot->FindBin(pitchX / 2))
				|| (iBin >= h_plot->FindBin(pitchX / 2 + epitchX)))
			h_innerPixel->SetBinContent(iBin, 0.0);
		else
			h_innerPixel->SetBinContent(iBin, h_innerPixel->GetMaximum());
	}
	h_innerPixel->SetFillColor(kBlue);
	h_innerPixel->SetLineColor(kGray);
	h_innerPixel->SetFillStyle(3004);
	h_innerPixel->SetLineWidth((Width_t) 0.5); // BJD: Really?!
	h_innerPixel->Draw("same");

	TF1* fit = tbutils::GetGausLineFit(histo, "turnon_fullrange", fitLeft,
			fitRight, 50.0, 9.0, 0.98);
	fit->SetLineColor(kRed);
	fit->SetLineWidth((Width_t) 2.3); // BJD: Really?!
	fit->Draw("L same");
	double xf90 = fit->GetX(0.9);
	double xf10 = fit->GetX(0.1);
	TPad *pad2 = new TPad("", "", 0.405, 0.27, 0.705, 0.52);
	pad2->Draw();
	pad2->cd();

	TPaveStats *paramBox = new TPaveStats(0.0, 0.0, 1.0, 1.0);
	paramBox->AddText("Active edge fit");
	char * line = new char[300];
	sprintf(line, "Pixel width=%.1f#pm%.1f #mum", pitchX - fit->GetParameter(0),
			fit->GetParError(0));
	paramBox->AddText(line);
	sprintf(line, "Resolution=%.1f#pm%.1f #mum", fit->GetParameter(1),
			fit->GetParError(1));
	paramBox->AddText(line);
	sprintf(line, "Plateau eff.=%.3f#pm%.3f", fit->GetParameter(2),
			fit->GetParError(2));
	paramBox->AddText(line);
	sprintf(line, "90%% (10%%)= %.1f (%.1f) #mum", pitchX - xf90,
			pitchX - xf10);
	paramBox->AddText(line);
	paramBox->SetFillColor(0);
	paramBox->SetTextColor(kBlack);
	paramBox->Draw();
	canvas->SaveAs(fileName);
	
	// added (2014-06-17) by Botho
	//config.saveToFile(this->name, histoName, (TNamed*) histo);
	config.saveToFile(this->name, histoName, canvas);

	delete canvas;
	delete fileName;
	delete fit;
	delete paramBox;
	return;

}

void EdgeEfficiency::ApplySensorMapStyle(const TbConfig& config, TH2D* histo,
		const char* title, const char* fileName, const char* histoName) const {
	TCanvas *canvas = new TCanvas("c1", "c1", 600, 504);
	histo->Draw("colz");
	histo->SetTitle(title);
	TStyle *tmpStyle = (TStyle*) gStyle->Clone("curstyle");
	gStyle->SetNdivisions(510, "x");
	gStyle->SetNdivisions(510, "y");
	gStyle->SetNdivisions(510, "z");
	gStyle->SetAxisColor(1, "x");
	gStyle->SetAxisColor(1, "y");
	gStyle->SetAxisColor(1, "z");
	gStyle->SetLabelColor(1, "x");
	gStyle->SetLabelColor(1, "y");
	gStyle->SetLabelColor(1, "z");
	gStyle->SetLabelFont(63, "x");
	gStyle->SetLabelFont(63, "y");
	gStyle->SetLabelFont(63, "z");
	gStyle->SetLabelOffset(0.007, "x");
	gStyle->SetLabelOffset(0.006, "y");
	gStyle->SetLabelOffset(0.005, "z");
	gStyle->SetLabelSize(16, "x");
	gStyle->SetLabelSize(16, "y");
	gStyle->SetLabelSize(16, "z");
	gStyle->SetTickLength(0.015, "x");
	gStyle->SetTickLength(0.015, "y");
	gStyle->SetTickLength(0.048, "z");
	gStyle->SetTitleOffset(0.8, "x");
	gStyle->SetTitleOffset(0.85, "y");
	gStyle->SetTitleOffset(1, "z");
	gStyle->SetTitleSize(22, "x");
	gStyle->SetTitleSize(22, "y");
	gStyle->SetTitleSize(16, "z");
	gStyle->SetTitleColor(1, "x");
	gStyle->SetTitleColor(1, "y");
	gStyle->SetTitleColor(1, "z");
	gStyle->SetTitleFont(53, "x");
	gStyle->SetTitleFont(53, "y");
	gStyle->SetTitleFont(63, "z");
	gStyle->SetBarWidth(1);
	gStyle->SetBarOffset(0);
	gStyle->SetDrawBorder(0);
	gStyle->SetOptLogx(0);
	gStyle->SetOptLogy(0);
	gStyle->SetOptLogz(0);
	gStyle->SetOptDate(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetOptTitle(kTRUE);
	gStyle->SetOptFit(0);
	gStyle->SetNumberContours(20);
	gStyle->GetAttDate()->SetTextFont(62);
	gStyle->GetAttDate()->SetTextSize(0.025);
	gStyle->GetAttDate()->SetTextAngle(0);
	gStyle->GetAttDate()->SetTextAlign(11);
	gStyle->GetAttDate()->SetTextColor(1);
	gStyle->SetDateX(0.01);
	gStyle->SetDateY(0.01);
	gStyle->SetEndErrorSize(2);
	gStyle->SetErrorX(0.5);
	gStyle->SetFuncColor(632);
	gStyle->SetFuncStyle(0);
	gStyle->SetFuncWidth(2);
	gStyle->SetGridColor(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetLegendBorderSize(4);
	gStyle->SetHatchesLineWidth(1);
	gStyle->SetHatchesSpacing(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameFillStyle(1);
	gStyle->SetFrameLineStyle(0);
	gStyle->SetFrameLineWidth(1);
	gStyle->SetFrameBorderSize(2);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetHistFillColor(0);
	gStyle->SetHistLineColor(1);
	gStyle->SetHistFillStyle(1001);
	gStyle->SetHistLineStyle(1);
	gStyle->SetHistLineWidth(1);
	gStyle->SetHistMinimumZero(kFALSE);
	gStyle->SetCanvasPreferGL(kFALSE);
	gStyle->SetCanvasColor(0);
	gStyle->SetCanvasBorderSize(2);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasDefH(500);
	gStyle->SetCanvasDefW(600);
	gStyle->SetCanvasDefX(10);
	gStyle->SetCanvasDefY(10);
	gStyle->SetPadColor(0);
	gStyle->SetPadBorderSize(2);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBottomMargin(0.10);
	gStyle->SetPadTopMargin(0.09);
	gStyle->SetPadLeftMargin(0.08);
	gStyle->SetPadRightMargin(0.11);
	gStyle->SetPadGridX(kFALSE);
	gStyle->SetPadGridY(kFALSE);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPaperSize(20, 26);
	gStyle->SetScreenFactor(1);
	gStyle->SetStatColor(0);
	gStyle->SetStatTextColor(1);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatFont(62);
	gStyle->SetStatFontSize(15);
	gStyle->SetStatStyle(1001);
	gStyle->SetStatFormat("6.4g");
	gStyle->SetStatX(1);
	gStyle->SetStatY(1);
	gStyle->SetStatW(0.45);
	gStyle->SetStatH(0.08);
	gStyle->SetStripDecimals(kTRUE);
	gStyle->SetTitleAlign(13);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFont(53);
	gStyle->SetTitleFontSize(0.048);
	gStyle->SetTitleStyle(1001);
	gStyle->SetTitleX(0.06);
	gStyle->SetTitleY(0.98);
	gStyle->SetTitleW(0);
	gStyle->SetTitleH(0);
	gStyle->SetLegoInnerR(0.5);
	gStyle->SetHeaderPS("");
	gStyle->SetTitlePS("");
	gStyle->SetFitFormat("5.4g");
	gStyle->SetPaintTextFormat("g");
	gStyle->SetLineScalePS(3);
	gStyle->SetColorModelPS(0);
	gStyle->SetTimeOffset(788918400);
	gStyle->SetLineColor(0);
	gStyle->SetLineStyle(1);
	gStyle->SetLineWidth(1);
	gStyle->SetFillColor(19);
	gStyle->SetFillStyle(1001);
	gStyle->SetMarkerColor(1);
	gStyle->SetMarkerSize(0.1);
	gStyle->SetMarkerStyle(20);
	gStyle->SetTextAlign(11);
	gStyle->SetTextAngle(0);
	gStyle->SetTextColor(1);
	gStyle->SetTextFont(62);
	gStyle->SetTextSize(0.05);

	histo->Draw("cont4");
	gStyle->SetOptTitle(kFALSE);
	canvas->UseCurrentStyle();
	canvas->Modified();
	canvas->Update();

	histo->Draw("colz");
	gPad->Update();
	TPaletteAxis *palette =
			(TPaletteAxis*) histo->GetListOfFunctions()->FindObject("palette");
	if (palette != 0) {
		palette->SetX1NDC(0.895);
		palette->SetX2NDC(0.935);
		palette->SetY1NDC(0.100);
		palette->SetY2NDC(0.910);
	}
	gStyle->SetOptTitle(kTRUE);
	canvas->UseCurrentStyle();
	canvas->Modified();
	canvas->Update();

	canvas->SaveAs(fileName);
	config.saveToFile(name, histoName, histo);
	gStyle = tmpStyle;
	delete canvas;
}

void EdgeEfficiency::ApplyPixelMapStyle(const TbConfig& config, TH2D* histo,
		const char* title, const char* fileName, const char* histoName) const {
	TCanvas *canvas = new TCanvas("c1", "c1", 1200, 300);
	histo->Draw("colz");
	histo->SetTitle(title);
	TStyle *tmpStyle = (TStyle*) gStyle->Clone("curstyle");
	gStyle->SetNdivisions(505, "x");
	gStyle->SetNdivisions(505, "y");
	gStyle->SetNdivisions(505, "z");
	gStyle->SetAxisColor(1, "x");
	gStyle->SetAxisColor(1, "y");
	gStyle->SetAxisColor(1, "z");
	gStyle->SetLabelColor(1, "x");
	gStyle->SetLabelColor(1, "y");
	gStyle->SetLabelColor(1, "z");
	gStyle->SetLabelFont(63, "x");
	gStyle->SetLabelFont(63, "y");
	gStyle->SetLabelFont(63, "z");
	gStyle->SetLabelOffset(0.007, "x");
	gStyle->SetLabelOffset(0.006, "y");
	gStyle->SetLabelOffset(0.005, "z");
	gStyle->SetLabelSize(16, "x");
	gStyle->SetLabelSize(16, "y");
	gStyle->SetLabelSize(16, "z");
	gStyle->SetTickLength(0.015, "x");
	gStyle->SetTickLength(0.015, "y");
	gStyle->SetTickLength(0.037, "z");
	gStyle->SetTitleOffset(1, "x");
	gStyle->SetTitleOffset(0.4, "y");
	gStyle->SetTitleOffset(1, "z");
	gStyle->SetTitleSize(0.058, "x");
	gStyle->SetTitleSize(16, "y");
	gStyle->SetTitleSize(16, "z");
	gStyle->SetTitleColor(1, "x");
	gStyle->SetTitleColor(1, "y");
	gStyle->SetTitleColor(1, "z");
	gStyle->SetTitleFont(53, "x");
	gStyle->SetTitleFont(53, "y");
	gStyle->SetTitleFont(63, "z");
	gStyle->SetBarWidth(1);
	gStyle->SetBarOffset(0);
	gStyle->SetDrawBorder(0);
	gStyle->SetOptLogx(0);
	gStyle->SetOptLogy(0);
	gStyle->SetOptLogz(0);
	gStyle->SetOptDate(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetOptTitle(kTRUE);
	gStyle->SetOptFit(1);
	gStyle->SetNumberContours(20);
	gStyle->GetAttDate()->SetTextFont(62);
	gStyle->GetAttDate()->SetTextSize(0.025);
	gStyle->GetAttDate()->SetTextAngle(0);
	gStyle->GetAttDate()->SetTextAlign(11);
	gStyle->GetAttDate()->SetTextColor(1);
	gStyle->SetDateX(0.01);
	gStyle->SetDateY(0.01);
	gStyle->SetEndErrorSize(2);
	gStyle->SetErrorX(0.5);
	gStyle->SetFuncColor(632);
	gStyle->SetFuncStyle(0);
	gStyle->SetFuncWidth(2);
	gStyle->SetGridColor(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetLegendBorderSize(4);
	gStyle->SetHatchesLineWidth(1);
	gStyle->SetHatchesSpacing(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameFillStyle(1);
	gStyle->SetFrameLineStyle(0);
	gStyle->SetFrameLineWidth(1);
	gStyle->SetFrameBorderSize(2);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetHistFillColor(0);
	gStyle->SetHistLineColor(1);
	gStyle->SetHistFillStyle(1001);
	gStyle->SetHistLineStyle(1);
	gStyle->SetHistLineWidth(1);
	gStyle->SetHistMinimumZero(kFALSE);
	gStyle->SetCanvasPreferGL(kFALSE);
	gStyle->SetCanvasColor(0);
	gStyle->SetCanvasBorderSize(2);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasDefH(500);
	gStyle->SetCanvasDefW(600);
	gStyle->SetCanvasDefX(10);
	gStyle->SetCanvasDefY(10);
	gStyle->SetPadColor(0);
	gStyle->SetPadBorderSize(2);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.12);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.12);
	gStyle->SetPadGridX(kFALSE);
	gStyle->SetPadGridY(kFALSE);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPaperSize(20, 26);
	gStyle->SetScreenFactor(1);
	gStyle->SetStatColor(0);
	gStyle->SetStatTextColor(1);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatFont(62);
	gStyle->SetStatFontSize(0.03);
	gStyle->SetStatStyle(1001);
	gStyle->SetStatFormat("6.4g");
	gStyle->SetStatX(0.995);
	gStyle->SetStatY(0.995);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.35);
	gStyle->SetStripDecimals(kTRUE);
	gStyle->SetTitleAlign(13);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFont(52);
	gStyle->SetTitleFontSize(0.08);
	gStyle->SetTitleStyle(1001);
	gStyle->SetTitleX(0.12);
	gStyle->SetTitleY(0.99);
	gStyle->SetTitleW(0);
	gStyle->SetTitleH(0);
	gStyle->SetLegoInnerR(0.5);
	gStyle->SetHeaderPS("");
	gStyle->SetTitlePS("");
	gStyle->SetFitFormat("5.4g");
	gStyle->SetPaintTextFormat("g");
	gStyle->SetLineScalePS(3);
	gStyle->SetColorModelPS(0);
	gStyle->SetTimeOffset(788918400);
	gStyle->SetLineColor(0);
	gStyle->SetLineStyle(1);
	gStyle->SetLineWidth(1);
	gStyle->SetFillColor(19);
	gStyle->SetFillStyle(1001);
	gStyle->SetMarkerColor(1);
	gStyle->SetMarkerSize(0.1);
	gStyle->SetMarkerStyle(20);
	gStyle->SetTextAlign(11);
	gStyle->SetTextAngle(0);
	gStyle->SetTextColor(1);
	gStyle->SetTextFont(62);
	gStyle->SetTextSize(0.05);

	histo->Draw("cont4");
	gStyle->SetOptTitle(kFALSE);
	canvas->UseCurrentStyle();
	canvas->Modified();
	canvas->Update();

	histo->Draw("colz");
	gPad->Update();
	TPaletteAxis *palette =
			(TPaletteAxis*) histo->GetListOfFunctions()->FindObject("palette");
	if (palette != 0) {
		palette->SetX1NDC(0.885);
		palette->SetX2NDC(0.913);
		palette->SetY1NDC(0.120);
		palette->SetY2NDC(0.880);
	}
	gStyle->SetOptTitle(kTRUE);
	canvas->UseCurrentStyle();
	canvas->Modified();
	canvas->Update();

	canvas->SaveAs(fileName);
	config.saveToFile(name, histoName, histo);
	gStyle = tmpStyle;
	delete canvas;
}
