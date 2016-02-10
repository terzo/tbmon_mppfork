#include "qshare2D.h"

void QShare2D::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();
	int maxToT = (int) floor(2304 / nCols); // gives reasonable range for both FE-I3 and FE-I4

	h_anyHitTrkMap = new TH2D("", ";Track X [#mum];Track Y [#mum]", pitchX / 2,
			0, pitchX * 2.0, pitchY / 2, 0, pitchY * 2.0);
	h_sharedHitTrkMap = new TH2D("", ";Track X [#mum];Track Y [#mum]",
			pitchX / 2, 0, pitchX * 2.0, pitchY / 2, 0, pitchY * 2.0);
	h_anyHitTrkMap_mirror = new TH2D("", ";Track X [#mum];Track Y [#mum]",
			pitchX / 2, 0, pitchX * 2.0, pitchY / 2, 0, pitchY * 2.0);
	h_sharedHitTrkMap_mirror = new TH2D("", ";Track X [#mum];Track Y [#mum]",
			pitchX / 2, 0, pitchX * 2.0, pitchY / 2, 0, pitchY * 2.0);
	h_directCharge = new TH3D("",
			";Track X [#mum];Track Y [#mum];Charge in Hit Pixel", pitchX / 2, 0,
			pitchX * 2.0, pitchY / 2, 0, pitchY * 2.0, maxToT, 0, maxToT);
	h_chargeShareScatter = new TH3D("",
			";Track X [#mum];Track Y [#mum];Fraction of Charge in Hit Pixel",
			pitchX / 2, 0, pitchX * 2.0, pitchY / 2, 0, pitchY * 2.0, 200, 0.0,
			1.0);

	numTracks = 0;
	doCharge = false;
	// range is set up for ToT... change if Q needed
	if (doCharge == true) {
		h_directCharge->GetZaxis()->SetRangeUser(0.0, 60000.0);
	}

}

void QShare2D::event(const TbConfig &config, const Event &event) {

	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		return;
	}
	// Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		return;
	}
	// Are we in the central region on the chip?
	if (event.fTrackCentralRegion != event::kGood) {
		return;
	}
	// Are we in bad regions or close to masked pixels?
	if (event.fTrackMaskedRegion != event::kGood) {
		return;
	}
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	// Check that we have only 1 cluster
	if (event.clusters.size() != 1) {
		return;
	}
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);
	if (matched == -1) {
		return;
	}
	// Require hits, but apparently not too many
	if (event.hits.size() == 0 or event.hits.size() > 4) {
		return;
	}
	numTracks++;

	// distance to a point 1/2 half a pixel from the edge and extending 
	// over two pixel lengths
	double pitchX = event.dut->getPitchX();
	double pitchY = event.dut->getPitchY();

	double trackXMod = fmod(event.trackX, 2 * pitchX);
	double trackYMod = fmod(event.trackY, 2 * pitchY);

	// Loop over all hits:
	bool hit_through = false;
	bool hit_vert = false;
	bool hit_hor = false;
	bool hit_diag = false;
	int hitCol = 0;
	int hitRow = 0;
	int hitToT = 0;
	double directCharge;

	for (int ii = 0; ii < event.hits.size(); ii++) {
		double diffX = fabs(event.hits.at(ii)->col * pitchX - event.trackX); // distance in x direction of track from centre of hit pixel
		double diffY = fabs(event.hits.at(ii)->row * pitchY - event.trackY); // distance in y direction of track from centre of hit pixel

		// Pixel firing has track through it
		if (diffX < 0.5 * pitchX and diffY < 0.5 * pitchY) {
			hit_through = true;
			hitCol = event.hits.at(ii)->col;
			hitRow = event.hits.at(ii)->row;
			hitToT = event.hits.at(ii)->tot;
		}
		// Pixel firing has track above or below
		if (diffX < 0.5 * pitchX and diffY > 0.5 * pitchY
				and diffY < 1.5 * pitchY) {
			hit_vert = true;
		}
		// Pixel firing has track to the left or right
		if (diffY < 0.5 * pitchY and diffX > 0.5 * pitchX
				and diffX < 1.5 * pitchX) {
			hit_hor = true;
		}
		// Pixel firing has track diagonaly to it
		if (diffX > 0.5 * pitchX and diffX < 1.5 * pitchX
				and diffY > 0.5 * pitchY and diffY < 1.5 * pitchY) {
			hit_diag = true;
		}
	}

	if (hit_through or hit_vert or hit_hor or hit_diag) {
		h_anyHitTrkMap->Fill(trackXMod, trackYMod);
		h_anyHitTrkMap_mirror->Fill(pitchX * 2.0 - trackXMod, trackYMod);
	}
	if (hit_through and (hit_vert or hit_hor or hit_diag)) {
		h_sharedHitTrkMap->Fill(trackXMod, trackYMod);
		h_sharedHitTrkMap_mirror->Fill(pitchX * 2.0 - trackXMod, trackYMod);
	}
	if (hit_through) {

		if (doCharge == true) {
			directCharge = event.dut->q(hitToT, hitCol, hitRow);
		} else {
			directCharge = hitToT;
		}
		h_directCharge->Fill(trackXMod, trackYMod, directCharge);

		if (doCharge == true) {
			if (hit_hor or hit_vert or hit_diag) {
				h_chargeShareScatter->Fill(trackXMod, trackYMod,
						directCharge
								/ cluster::getSumChargePP(
										event.clusters.at(matched), event));
			}
		} else {
			if (hit_hor or hit_vert or hit_diag) {
				h_chargeShareScatter->Fill(trackXMod, trackYMod,
						directCharge
								/ cluster::getSumTot(
										event.clusters.at(matched)));
			}
		}

	} else {
		h_chargeShareScatter->Fill(trackXMod, trackYMod, 0.0);
	}

}

void QShare2D::finalize(const TbConfig &config) {

	// Get overall charge sharing probability
	double nAllHits = h_anyHitTrkMap->GetEntries();
	double nSharedHits = h_sharedHitTrkMap->GetEntries();

	TBALOG(kINFO) << "Tracks:                 " << numTracks << endl;
	TBALOG(kINFO) << "All Hits:               " << nAllHits << endl;
	TBALOG(kINFO) << "Shared Hits:            " << nSharedHits << endl;
	TBALOG(kINFO) << "Overall Charge Sharing: " << nSharedHits / nAllHits
			<< endl;

	// not exactly sure what this is for... but not touching it (BJD)
	for (int i = 0; i < h_anyHitTrkMap->GetNbinsX(); i++) {
		for (int j = 0; j < h_anyHitTrkMap->GetNbinsY(); j++) {
			if (h_anyHitTrkMap->GetBinContent(i, j) == 0) {
				h_anyHitTrkMap->SetBinContent(i, j, 1000000.0);
			}
		}
	}

	config.drawAndSave(name, "anyHitTrkMap", "colz", h_anyHitTrkMap);
	config.drawAndSave(name, "sharedHitTrkMap", "colz", h_sharedHitTrkMap);
	config.drawAndSave(name, "anyHitTrkMap_mirror", "colz",
			h_anyHitTrkMap_mirror);
	config.drawAndSave(name, "sharedHitTrkMap_mirror", "colz",
			h_sharedHitTrkMap_mirror);
	config.drawAndSave(name, "directCharge", h_directCharge);
	config.drawAndSave(name, "chargeShareScatter", h_chargeShareScatter);
	config.drawAndSave(name, "directChargeProfile",
			h_directCharge->Project3DProfile());
	config.drawAndSave(name, "sharedChargeProfile",
			h_chargeShareScatter->Project3DProfile());

	drawAspect(config, name, "ChargeSharingProb2DAspect", h_anyHitTrkMap,
			h_sharedHitTrkMap);

	h_sharedHitTrkMap->Divide(h_anyHitTrkMap);
	h_sharedHitTrkMap->GetZaxis()->SetRangeUser(0, 1);
	h_sharedHitTrkMap->SetStats(kFALSE);
	config.drawAndSave(name, "ChargeSharing2D", "cont4z", h_sharedHitTrkMap);

	h_sharedHitTrkMap_mirror->Divide(h_anyHitTrkMap_mirror);
	h_sharedHitTrkMap_mirror->GetZaxis()->SetRangeUser(0, 1);
	h_sharedHitTrkMap_mirror->SetStats(kFALSE);
	config.drawAndSave(name, "ChargeSharing2D_mirror", "cont4z",
			h_sharedHitTrkMap_mirror);

}

void QShare2D::drawAspect(const TbConfig& config, const char* analysisName,
		const char* histoName, TH2D* h_any, TH2D* h_shared) const {

	char* fileName = config.buildHistName(analysisName, histoName);
	const DUT* dut = config.getDut(this->iden);

	// Get pixel efficiency map  
	TH2D* hr = (TH2D*) h_shared->Clone(TString::Format("hr"));
	hr->SetTitle("");
	hr->SetStats(false);
	hr->Divide(h_shared, h_any, 1.0, 1.0, "B");

	TH2D* hrA = (TH2D*) hr->Clone(TString::Format("hrA"));

	// In order make the palette readable I make two plots with different aspect ratios
	// one which I use only to get the palette
	// use a fixed palette in order to compare with other sensors

	//Remove the ticks on the axis
	//hrA->GetXaxis()->SetTickLength(0.00000000001);
	hrA->GetYaxis()->SetTickLength(0.01);
	//Reduce the axis label entries
	hrA->GetXaxis()->SetNdivisions(4);
	hrA->GetYaxis()->SetNdivisions(2);
	hrA->GetXaxis()->SetLabelOffset(0.005);
	hrA->GetYaxis()->SetLabelOffset(0.01);
	hrA->GetXaxis()->SetLabelSize(0.12);
	hrA->GetYaxis()->SetLabelSize(0.12);
	hrA->SetTitle(";;");

	// Draw the palette separately, since the aspect ratio will mess it up otherwise
	TCanvas* c_temp = new TCanvas("c_temp", "", 1200, 1200);
	c_temp->cd();
	gStyle->SetPalette(1, 0);
	hrA->Draw("cont4z");
	hrA->GetZaxis()->SetRangeUser(0, 1);
	gPad->Update();
	TPaletteAxis* palette =
			(TPaletteAxis*) hrA->GetListOfFunctions()->FindObject("palette");
	if (palette != 0) {
		palette->SetX1NDC(0.4);
		palette->SetX2NDC(0.6);
		palette->SetY1NDC(0.03);
		palette->SetY2NDC(0.97);

		TCanvas* c_palette = new TCanvas("c_palette", "", 1200, 1200);
		c_palette->cd();
		palette->Draw();
		char* histoNamePalette = new char[600];
		sprintf(histoNamePalette, "%s-palette", histoName);
		char* fileNamePalette = config.buildHistName(analysisName,
				histoNamePalette);
		c_palette->SaveAs(fileNamePalette);

		delete c_palette;
		delete histoNamePalette;
		delete fileNamePalette;
	} else {
		TBALOG(kERROR) << "Palette not found!" << endl;
	}
	delete c_temp;

	//This is canvas with correct aspect ratio
	double k = dut->getPitchY() / dut->getPitchX();
	TCanvas* canvas = new TCanvas("c1", "c1", 1200, floor(1200 * k));
	canvas->cd();
	gStyle->SetPalette(1, 0);
	hrA->GetZaxis()->SetRangeUser(0, 1);
	hrA->Draw("cont4");
	gPad->Update();
	canvas->SaveAs(fileName);
	//config.saveToFile(name,histoName,hrA);

	delete canvas;
	delete fileName;
	delete hrA;

}
