#include "edgeefficiencyshift.h"

void EdgeEfficiencyShift::init(TbConfig &config) {
	const DUT* dut = config.getDut(this->iden);
	ntracks = 0;
	ntracksedge = 0;
	nhitsedge = 0;

	leftTrackPixelMapEdge = new TH2D("",
			"Left Track Pixel Map Edge;Long side [#mum];Short Side[#mum]", 100,
			-2 * dut->getPitchX(), dut->getPitchX(), 2500, 0,
			dut->getPitchY() * 165);
	leftHitPixelMapEdge = new TH2D("",
			"Left Hit Pixel Map Edge;Long side [#mum];Short Side[#mum]", 100,
			-2 * dut->getPitchX(), dut->getPitchX(), 2500, 0,
			dut->getPitchY() * 165);
	rightTrackPixelMapEdge = new TH2D("",
			"right Track Pixel Map Edge;Long side [#mum];Short Side[#mum]", 100,
			6400, 7600, 2500, 0, dut->getPitchY() * 165);
	rightHitPixelMapEdge = new TH2D("",
			"right Hit Pixel Map Edge;Long side [#mum];Short Side[#mum]", 100,
			6400, 7600, 2500, 0, dut->getPitchY() * 165);
	rightHitOverlay = new TH2D("",
			"right Hit Overlay;Long side [#mum];Short Side[#mum]", 100, 6400,
			7600, 16, 0, dut->getPitchY() * 16);
	rightTrackOverlay = new TH2D("",
			"right Track Overlay;Long side [#mum];Short Side[#mum]", 100, 6400,
			7600, 16, 0, dut->getPitchY() * 16);
	hitOverlay = new TH2D("", "Hit Overlay;Long side [#mum];Short Side[#mum]",
			100, 0, 1200, 16, 0, dut->getPitchY() * 16);
	trackOverlay = new TH2D("",
			"Track Overlay;Long side [#mum];Short Side[#mum]", 100, 0, 1200, 16,
			0, dut->getPitchY() * 16);

	//charge stuff "borrowed" from qEfficiency
	h_chargeScatterRight = new TH3D("", "", 100, 6400, 7600, 16, 0,
			16 * dut->getPitchY(), 120, 0, 60e3);
	h_chargeScatterFolded = new TH3D("", "", 100, 0, 1200, 16, 0,
			16 * dut->getPitchY(), 120, 0, 60e3);
}

void EdgeEfficiencyShift::event(const TbConfig &config, const Event &event) {
	//Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		return;
	}
	//Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		return;
	}

	ntracks++;

	TBALOG( kDEBUG2 ) << "Event: " << config.currentEntry << " Run: "
			<< config.currentRun << endl;
	TBALOG( kDEBUG2 ) << "Good track." << endl;

	TBALOG( kDEBUG2 ) << "Edge track, col="
			<< tbutils::getCol(event.trackX, event) << " row="
			<< tbutils::getRow(event.trackY, event) << endl;

	//Make sure we are not too far off the sensor in the long (X) direction for edge studies

	if (event.trackX < -2.0 * event.dut->getPitchX() || event.trackX > 7600) {
		return;
	}
	if (event.trackX > 1.0 * event.dut->getPitchX() && event.trackX < 6400) {
		return;
	}

	//The edge in Y (ganged- and inter-ganged pixels) should still be ignored
	//if( (event.trackY + 0.5*event.dut->getPitchY()) < (event.dut->getPitchY() * event::skipRows) ) {return;}
	//if( (event.trackY + 0.5*event.dut->getPitchY()) > (event.dut->getPitchY() * (event::nrows - event::skipRows)) ) {return;}

	//Count the tracks
	ntracksedge++;

	TBALOG( kDEBUG2 ) << "Event: " << config.currentEntry << " Run: "
			<< config.currentRun << endl;
	TBALOG( kDEBUG2 ) << "Accepted track." << endl;

	//overlay the middle 8 pixels of a shifted set of 10 pixels
	int trackrow = tbutils::getRow(event.trackY, event);
	bool edgeofset = false;
	if (trackrow % 10 == 0 || trackrow % 10 == 9)
		edgeofset = true;
	double overlayY = (trackrow / 10) * event.dut->getPitchY()
			+ fmod(event.trackY + event.dut->getPitchY() * 0.5,
					event.dut->getPitchY());

	//3400um is the middle of the chip in x direction
	double foldedEdgeX =
			(event.trackX < 3400) ? (6800 - event.trackX) : event.trackX;
	foldedEdgeX -= 6400;

	{

		leftTrackPixelMapEdge->Fill(event.trackX, event.trackY);
		rightTrackPixelMapEdge->Fill(event.trackX, event.trackY);

		if (!edgeofset) {
			rightTrackOverlay->Fill(event.trackX, overlayY);
			trackOverlay->Fill(foldedEdgeX, overlayY);
		}
	}

	//Are there any hits in the DUT
	if (event.hits.size() > 0) {
		TBALOG( kDEBUG2 ) << "Found hit in DUT!" << endl;
	}

	//Look for matched clusters
	bool cluMatch(false);
	int clusterIndex(-1);
	if (event.trackX > 200 && event.trackX < 6600) {
		clusterIndex = cluster::getMatched(event.clusters, event);
		if (clusterIndex == -1) {
			return;
		}
	} else {
		clusterIndex = getEdgeMatched(event.clusters, event);
		if (clusterIndex == -1) {
			return;
		}
	}
	TBALOG( kDEBUG2 ) << "Fount matching cluster!" << endl;

	nhitsedge++;
	{

		leftHitPixelMapEdge->Fill(event.trackX, event.trackY);
		rightHitPixelMapEdge->Fill(event.trackX, event.trackY);

		if (!edgeofset) {
			rightHitOverlay->Fill(event.trackX, overlayY);
			hitOverlay->Fill(foldedEdgeX, overlayY);
		}

		//charge
		if (!edgeofset) {
			//h_chargeScatterRight->Fill( event.trackX,overlayY, cluster::getSumChargePP(event.clusters[clusterIndex],event));
			//h_chargeScatterFolded->Fill( foldedEdgeX,overlayY, cluster::getSumChargePP(event.clusters[clusterIndex],event));
			h_chargeScatterRight->Fill(event.trackX, overlayY,
					cluster::getSumTot(event.clusters[clusterIndex]));
			h_chargeScatterFolded->Fill(foldedEdgeX, overlayY,
					cluster::getSumTot(event.clusters[clusterIndex]));
		}
	}

}

void EdgeEfficiencyShift::finalize(const TbConfig &config) {
	const DUT* dut = config.getDut(this->iden);
	//TStyle* style = (TStyle*)gStyle->Clone("curstyle");
	tbutils::HistStyle();

	//config.drawToFile(name, "tracks", "cont4", trackMap);
	//config.drawToFile(name, "tracksedge", "colz", trackMapEdge);
	//config.drawToFile(name, "hitsanyedge", "colz", anyhitMapEdge);
	//config.drawToFile(name, "hitsedge", "colz", hitMapEdge);
	//config.drawToFile(name, "lefttrackspixeledge", "colz", leftTrackPixelMapEdge);
	//config.drawToFile(name, "hitspixeledge", "colz", leftHitPixelMapEdge);

	config.saveToFile(name, "lefttrackspixeledge", leftTrackPixelMapEdge);
	config.saveToFile(name, "lefthitspixeledge", leftHitPixelMapEdge);
	config.saveToFile(name, "righttrackspixeledge", rightTrackPixelMapEdge);
	config.saveToFile(name, "righthitspixeledge", rightHitPixelMapEdge);
	config.saveToFile(name, "righthitoverlay", rightHitOverlay);
	config.saveToFile(name, "righttrackoverlay", rightTrackOverlay);
	config.saveToFile(name, "hitoverlay", hitOverlay);
	config.saveToFile(name, "trackoverlay", trackOverlay);
	TBALOG( kINFO ) << " Tracks:             " << ntracks << endl;
	TBALOG( kINFO ) << " Edge tracks:        " << ntracksedge << endl;
	TBALOG( kINFO ) << " Edge tracks w/ hit: " << nhitsedge << endl;

	//charge stuff
	config.saveToFile(this->name, "chargeScatterRight", h_chargeScatterRight);
	TProfile2D* p_chargeProfileRight = h_chargeScatterRight->Project3DProfile(
			"yx");
	config.saveToFile(this->name, "chargeProfileRight", p_chargeProfileRight);

	config.saveToFile(this->name, "chargeScatterFolded", h_chargeScatterFolded);
	TProfile2D* p_chargeProfileFolded = h_chargeScatterFolded->Project3DProfile(
			"yx");
	config.saveToFile(this->name, "chargeProfileFolded", p_chargeProfileFolded);

	//get profiles for individual sets of shifted pixels (kindly stolen from qEfficiency.cc)
	char sliceName[20], landauSliceName[20];
	TH2D* preProjectSlice;
	TProfile* p_Slice;
	TH3D* chargetemp;
	for (int i = 0; i < 15; i++) {
		h_chargeScatterFolded->GetYaxis()->SetRangeUser(dut->getPitchY() * i,
				dut->getPitchY() * (i + 1));
		preProjectSlice = (TH2D*) h_chargeScatterFolded->Project3D("zx");
		sprintf(sliceName, "chargeSliceSet%u", i);
		sprintf(landauSliceName, "landauSliceSet%u", i);
		p_Slice = preProjectSlice->ProfileX("", 1, -1, "");
		config.saveToFile(this->name, sliceName, p_Slice);
		delete preProjectSlice;
		delete p_Slice;
		//clone the charge scatter
		chargetemp = new TH3D(*h_chargeScatterFolded);
		config.saveToFile(this->name, landauSliceName,
				makeLandauProfile(config, chargetemp));
		delete chargetemp;
	}

	//Reset axis
	h_chargeScatterFolded->GetYaxis()->SetRange(1,
			h_chargeScatterFolded->GetYaxis()->GetNbins());

	//efficiency maps  
	TH2D* leftEffPixelMapEdge = (TH2D*) leftHitPixelMapEdge->Clone(
			TString::Format("leftEffPixelMapEdge"));
	leftEffPixelMapEdge->SetTitle("");
	leftEffPixelMapEdge->SetStats(false);
	leftEffPixelMapEdge->Divide(leftHitPixelMapEdge, leftTrackPixelMapEdge, 1.,
			1., "B");
	config.saveToFile(name, "lefteffpixelmapedge", leftEffPixelMapEdge);
	drawToFileSensorAspect(config, name, "lefteffpixelmapedgeaspect", "colz",
			leftEffPixelMapEdge);

	TH2D* rightEffPixelMapEdge = (TH2D*) rightHitPixelMapEdge->Clone(
			TString::Format("rightEffPixelMapEdge"));
	rightEffPixelMapEdge->SetTitle("");
	rightEffPixelMapEdge->SetStats(false);
	rightEffPixelMapEdge->Divide(rightHitPixelMapEdge, rightTrackPixelMapEdge,
			1., 1., "B");
	config.saveToFile(name, "righteffpixelmapedge", rightEffPixelMapEdge);
	drawToFileSensorAspect(config, name, "righteffpixelmapedgeaspect", "colz",
			rightEffPixelMapEdge);

	TH2D* rightEffOverlay = (TH2D*) rightHitOverlay->Clone(
			TString::Format("rightEffOverlay"));
	rightEffOverlay->SetTitle("");
	rightEffOverlay->SetStats(false);
	rightEffOverlay->Divide(rightHitOverlay, rightTrackOverlay, 1., 1., "B");
	config.saveToFile(name, "righteffoverlay", rightEffOverlay);
	drawToFileSensorAspect(config, name, "righteffoverlayaspect", "colz",
			rightEffOverlay);

	TH2D* effOverlay = (TH2D*) hitOverlay->Clone(TString::Format("effOverlay"));
	effOverlay->SetTitle("");
	effOverlay->SetStats(false);
	effOverlay->Divide(hitOverlay, trackOverlay, 1., 1., "B");
	config.saveToFile(name, "effoverlay", effOverlay);
	drawToFileSensorAspect(config, name, "effoverlayaspect", "colz",
			effOverlay);

	//Reset the ROOT style
	gROOT->SetStyle();

}

void EdgeEfficiencyShift::drawToFileSensorAspect(const TbConfig& config,
		const char* analysisName, const char* histoName, const char* drawOpts,
		TH1* h1, TH1* h2, TH1* h3, TH1* h4) const {
	const DUT* dut = config.getDut(this->iden);
	//gets the correct aspect ratio of the pixel
	char* fileName = config.buildHistName(analysisName, histoName);
	double k = dut->getPitchY() / dut->getPitchX();
	TCanvas* canvas = new TCanvas("", "", 1200, floor(1200 * k));
	canvas->cd();
	h1->Draw(drawOpts);
	if (h2 != NULL) {
		h2->Draw("same");
	}
	if (h3 != NULL) {
		h3->Draw("same");
	}
	if (h4 != NULL) {
		h4->Draw("same");
	}
	canvas->SaveAs(fileName);
	delete canvas;
	delete[] fileName;
}

void EdgeEfficiencyShift::drawToFileEdgeEff(const TbConfig& config,
		const Event &event, const char* analysisName, const char* histoName,
		TObject* histo) const {
	char* fileName = config.buildHistName(analysisName, histoName);
	//make sure both a graph and hist can be used in the fit
	TGraphAsymmErrors* grhisto = 0;
	TH1D* hhisto = 0;
	bool isGraph = false;
	if (histo->InheritsFrom("TGraph")) {
		isGraph = true;
		grhisto = (TGraphAsymmErrors*) histo;
	} else if (histo->InheritsFrom("TH1D")) {
		hhisto = (TH1D*) histo;
	} else {
		TBALOG( kERROR )
				<< "You have to have either a TH1D or TGraph for the edge fit"
				<< endl;
		exit(1);
	}
	assert(tbutils::XOR(grhisto, hhisto));
	TCanvas* canvas = new TCanvas("", "", 1200, 300);
	canvas->cd();
	TString title = ";Long pixel[ #mum ];Active edge efficiency";
	TH1* hplot = 0;
	if (isGraph) {
		grhisto->SetLineWidth((Width_t) 1.3);
		grhisto->SetMarkerStyle(20);
		grhisto->SetMarkerSize(0.7);
		grhisto->Draw("AP");
		grhisto->SetTitle(title);
		grhisto->GetXaxis()->SetLabelSize(0.07);
		grhisto->GetYaxis()->SetLabelSize(0.07);
		grhisto->GetYaxis()->SetTitleSize(0.08);
		grhisto->GetXaxis()->SetTitleSize(0.08);
		grhisto->GetYaxis()->SetTitleOffset(0.3);
		grhisto->GetXaxis()->SetTitleOffset(0.9);
		//grhisto->GetXaxis()->SetTitleSize(0.09);
		hplot = grhisto->GetHistogram();
		hplot->SetFillColor(kOrange);
		hplot->Draw("same");
	} else {
		hhisto->SetLineWidth((Width_t) 1.3);
		hhisto->SetTitle(title);
		grhisto->GetYaxis()->SetTitleSize(0.09);
		grhisto->GetYaxis()->SetTitleOffset(0.03);
		//grhisto->GetXaxis()->SetTitleSize(0.09);
		hhisto->Draw();
		hplot = hhisto;
	}

	//Draw an area covering the next to last pixel for visualization
	int outeredgebin = hplot->FindBin(600.0);
	TH1D* h_innerpixel = (TH1D*) hplot->Clone("h_innerpixel");
	for (int ibin = 1; ibin != h_innerpixel->GetNbinsX() + 1; ++ibin) {
		if (ibin < outeredgebin)
			h_innerpixel->SetBinContent(ibin, 0.);
		else
			h_innerpixel->SetBinContent(ibin, h_innerpixel->GetMaximum());
	}
	h_innerpixel->SetFillColor(kBlue);
	h_innerpixel->SetFillStyle(3004);
	h_innerpixel->SetLineWidth(0.5);
	//h_innerpixel->SetLineColor(kGray);
	h_innerpixel->Draw("same");

	tbutils::myMultiLineText(0.7, 0.26, "Inner|pixel", kWhite, kBlue + 2);
	tbutils::myMultiLineText(0.64, 0.26, "Edge|pixel", kWhite, kBlack);
	//Redraw the plot to get it above
//   if(isGraph) {
//     grhisto->Draw("ALP,same");
//   } else {
//     hhisto->Draw("same");
//   }

	TF1* fit = tbutils::GetGausLineFit(histo, "turnon_fullrange", 10, 600, 50.,
			9.0, 0.98);
	fit->SetLineColor(kRed);
	fit->SetLineWidth((Width_t) 2.3);
	fit->Draw("L,same");
	double xf90 = fit->GetX(0.9);
	double xf10 = fit->GetX(0.1);
	TPad *pad2 = new TPad("", "", 0.28, 0.27, 0.5, 0.51);
	pad2->Draw();
	pad2->cd();
	//myText();

	TPaveStats *paramBox = new TPaveStats(0.0, 0.0, 1.0, 1.0);
	//paramBox->AddText("Active edge fit");
	char * line = new char[300];
	sprintf(line, "Pixel width=%.1f#pm%.1f #mum", 600.0 - fit->GetParameter(0),
			fit->GetParError(0));
	paramBox->AddText(line);
	sprintf(line, "Resolution=%.1f#pm%.1f #mum", fit->GetParameter(1),
			fit->GetParError(1));
	paramBox->AddText(line);
	sprintf(line, "Plateau eff.=%.3f#pm%.3f", fit->GetParameter(2),
			fit->GetParError(2));
	paramBox->AddText(line);
	sprintf(line, "90%% (10%%)= %.1f(%.1f) #mum", 600.0 - xf90, 600.0 - xf10);
	paramBox->AddText(line);
	paramBox->SetFillColor(0);
	paramBox->SetTextColor(kBlack);
	paramBox->Draw();
	canvas->SaveAs(fileName);
	TString rootfilename(fileName);
	rootfilename.ReplaceAll(".", "_");
	rootfilename.ReplaceAll("-", "_");
	rootfilename.ReplaceAll("/", "_");
	canvas->SetName(rootfilename.Data());
	rootfilename = fileName;
	rootfilename.ReplaceAll(".", "_");
	rootfilename += ".C";
	canvas->SaveAs(rootfilename.Data());
	delete canvas;
	delete fileName;
	delete fit;
	delete paramBox;
	//delete pad2;

	return;
}

void EdgeEfficiencyShift::drawToFilePixelMapEff(const TbConfig& config,
		const Event &event, const char* analysisName, const char* histoName,
		TH2D* hitPixelMap, TH2D* trackPixelMap) const {
	char* fileName = config.buildHistName(analysisName, histoName);
	trackPixelMap->Rebin2D(8, 2);
	hitPixelMap->Rebin2D(8, 2);
	//Find the pixel efficiency map  
	TH2D* effPixelMap = (TH2D*) hitPixelMap->Clone(
			TString::Format("effPixelMap"));
	effPixelMap->SetTitle("");
	effPixelMap->SetStats(false);
	effPixelMap->Divide(hitPixelMap, trackPixelMap, 1., 1., "B");

	TH2D* effPixelMapAspect = (TH2D*) effPixelMap->Clone(
			TString::Format("effPixelMapAspect"));

	// In order make the palette readable I make two plots with different aspect ratios
	// one which I use only to get the palette
	// use a fixed palette in order to compare with other sensors

	//Remove the ticks on the axis
	//effPixelMapAspect->GetXaxis()->SetTickLength(0.00000000001);
	effPixelMapAspect->GetYaxis()->SetTickLength(0.01);
	//Reduce the axis label entries
	effPixelMapAspect->GetXaxis()->SetNdivisions(4);
	effPixelMapAspect->GetYaxis()->SetNdivisions(2);
	effPixelMapAspect->GetXaxis()->SetLabelOffset(0.045);
	effPixelMapAspect->GetYaxis()->SetLabelOffset(0.012);
	effPixelMapAspect->GetXaxis()->SetLabelSize(0.12);
	effPixelMapAspect->GetYaxis()->SetLabelSize(0.12);
	//Remove the titles
	effPixelMapAspect->SetTitle(";;");

	//This is canvas with correct aspect ratio
	double k = event.dut->getPitchY() / event.dut->getPitchX();
	TCanvas* canvas = new TCanvas("", "", 1600, floor(1600 * k));
	canvas->cd();
	gStyle->SetPalette(1, 0);
	effPixelMapAspect->Draw("colz");
	gPad->Update();
	TPaletteAxis *palette =
			(TPaletteAxis*) effPixelMapAspect->GetListOfFunctions()->FindObject(
					"palette");
	//Remove the Palette
	if (palette != 0) {
		palette->SetY1NDC(1.0);
		palette->SetY2NDC(0.9999999);
		palette->SetX1NDC(1.0);
		palette->SetX2NDC(0.9999999);
	} else {
		TBALOG( kERROR ) << "palette not found!" << endl;
		//exit(-1);    
	}
	canvas->SaveAs(fileName);
	config.saveToFile(name, histoName, effPixelMapAspect);

	delete canvas;
	delete fileName;
	delete effPixelMapAspect;

	//This is canvas with wrong aspect ratio used to get the large palette
	TCanvas* canvas2 = new TCanvas("", "", 1600, 1600);
	canvas2->cd();
	gStyle->SetPalette(1, 0);
	effPixelMap->Draw("colz");
	gPad->Update();
	palette = (TPaletteAxis*) effPixelMap->GetListOfFunctions()->FindObject(
			"palette");
	if (palette != 0) {
		palette->SetY2NDC(0.99);
		palette->SetY1NDC(0.03);
		palette->SetX2NDC(0.7);
	} else {
		TBALOG( kERROR ) << "palette not found!" << endl;
		//exit(-1);    
	}
	char* histoName2 = new char[600];
	sprintf(histoName2, "%s-large", histoName);
	char* fileName2 = config.buildHistName(analysisName, histoName2);
	canvas2->SaveAs(fileName2);
	delete fileName2;
	//delete canvas2;

	if (palette != 0) {
		TCanvas* canvas3 = new TCanvas("", "", 1600, 1600);
		canvas3->cd();
		palette->Draw();

		char* histoName3 = new char[600];
		sprintf(histoName3, "%s-largepalette", histoName);
		char* fileName3 = config.buildHistName(analysisName, histoName3);
		canvas3->SaveAs(fileName3);
		delete canvas3;
	}
	//delete canvas2;
	delete effPixelMap;

}

bool EdgeEfficiencyShift::isEdgeMatched(const event::cluster_t &clu,
		const Event& event) {
	bool matched(false);
	int trackrow = tbutils::getRow(event.trackY, event);
	//int trackcol = tbutils::getCol(event.trackX);

	//if outside of chip, use only col 0
	if (event.trackX < -400) {
		for (vector<PllHit*>::const_iterator i = clu.begin(); i != clu.end();
				i++) {
			if (((*i)->col == 0)
					&& ((*i)->row < trackrow + 2 && (*i)->row > trackrow - 2)) {
				matched = true;
				break;
			}
		}
	}
	//if in col 0, use col 0 and col 1
	else if (event.trackX < 200) {
		for (vector<PllHit*>::const_iterator i = clu.begin(); i != clu.end();
				i++) {
			if (((*i)->col < 2)
					&& ((*i)->row < trackrow + 2 && (*i)->row > trackrow - 2)) {
				matched = true;
				break;
			}
		}
	}
	//if outside of chip, use only col 17
	else if (event.trackX > 7200) {
		for (vector<PllHit*>::const_iterator i = clu.begin(); i != clu.end();
				i++) {
			if (((*i)->col == 17)
					&& ((*i)->row < trackrow + 2 && (*i)->row > trackrow - 2)) {
				matched = true;
				break;
			}
		}
	}
	//if in col 17, use col 16 and col 17
	else if (event.trackX > 6600) {
		for (vector<PllHit*>::const_iterator i = clu.begin(); i != clu.end();
				i++) {
			if (((*i)->col > 15)
					&& ((*i)->row < trackrow + 2 && (*i)->row > trackrow - 2)) {
				matched = true;
				break;
			}
		}
	}
	return matched;
}

int EdgeEfficiencyShift::getEdgeMatched(
		const vector<event::cluster_t> &clusters, const Event &event) {
	int index(-1);
	for (int ii = 0; ii < clusters.size(); ii++) {
		if (isEdgeMatched(clusters.at(ii), event)) {
			index = ii;
			break;
		}
	}
	return (index);
}

//again with many thanks to the skilled authors of qEfficiency!
TProfile* EdgeEfficiencyShift::makeLandauProfile(const TbConfig& config,
		TH3D* h_scatter, int rebin) {
	TH2D* h_slice = NULL;
	TH1D* h_sliceLandau = NULL;

	int nBins = ceil(h_scatter->GetXaxis()->GetNbins() / (double) rebin);

	TProfile* p_landauProfile = new TProfile("", "", nBins,
			h_scatter->GetXaxis()->GetXmin(), h_scatter->GetXaxis()->GetXmax());
	double* landauConst = new double[nBins];
	double* landauMPV = new double[nBins];
	double* landauSigma = new double[nBins];

	for (int i = 1; i < nBins; i++) {
		h_scatter->GetXaxis()->SetRange(i * rebin,
				(i + 1) * rebin <= h_scatter->GetXaxis()->GetNbins() ?
						(i + 1) * rebin : h_scatter->GetXaxis()->GetNbins());

		h_slice = (TH2D*) h_scatter->Project3D("zy");
		h_sliceLandau = h_slice->ProjectionY();
		TFitResultPtr fitRes = h_sliceLandau->Fit("landau", "SQ");
		if (!fitRes->IsValid()) {
			delete h_sliceLandau;
			delete h_slice;
			continue;
		}
		landauConst[i] = fitRes->Value(0);
		landauMPV[i] = fitRes->Value(1);
		landauSigma[i] = fitRes->Value(2);

		if (config.logLevel >= kDEBUG) {
			cout << "Fit result [i] = " << landauConst[i] << ", "
					<< landauMPV[i] << ", " << landauSigma[i] << endl;
		}

		p_landauProfile->SetBinContent(i, landauMPV[i]);
		p_landauProfile->SetBinError(i, landauSigma[i]);
		p_landauProfile->SetBinEntries(i, h_sliceLandau->GetEntries());

		delete h_sliceLandau;
		delete h_slice;
	}

	//reset x axis
	h_scatter->GetXaxis()->SetRange(1, h_scatter->GetXaxis()->GetNbins());

	return p_landauProfile;
}

//cleans up stray one entry bins in the charge profile
TProfile2D* EdgeEfficiencyShift::cleanUp(const TH3D& chargeScatter) {
	int nent;
	int ent_limit = 1;

	TProfile2D* p_chargeProfile = chargeScatter.Project3DProfile("yx");

	for (int i = 1; i <= chargeScatter.GetXaxis()->GetNbins(); i++) {
		for (int j = 1; j <= chargeScatter.GetYaxis()->GetNbins(); j++) {
			nent = 0;
			for (int k = 1; k <= chargeScatter.GetZaxis()->GetNbins(); k++) {
				nent += (int) chargeScatter.GetBinContent(i, j, k);
			}
			if (nent <= ent_limit)
				p_chargeProfile->SetBinContent(i, j, 0.0);
		}
	}

	return p_chargeProfile;
}
;

