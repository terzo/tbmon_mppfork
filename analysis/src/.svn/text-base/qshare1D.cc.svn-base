#include "qshare1D.h"

void QShare1D::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();
	// hack, but gives reasonable range for both FE-I3 and FE-I4
	int maxToT = (int) floor(2304 / nCols);

	h_tracks = new TH1D("", ";Y-Position in Pixel [#mum];Tracks",
			(int) (0.5 * pitchY), 0, pitchY);
	h_allHitsCenter = new TH1D("", ";Y-Position in Pixel [#mum]",
			(int) (0.5 * pitchY), 0, pitchY);
	h_nHits1Center = new TH1D("", ";Y-Position in Pixel [#mum]",
			(int) (0.5 * pitchY), 0, pitchY);
	h_nHits2Center = new TH1D("", ";Y-Position in Pixel [#mum]",
			(int) (0.5 * pitchY), 0, pitchY);
	h_nHits3PlusCenter = new TH1D("", ";Y-Position in Pixel [#mum]",
			(int) (0.5 * pitchY), 0, pitchY);
	h_allHitsNotCenter = new TH1D("", ";Y-Position in Pixel [#mum]",
			(int) (0.5 * pitchY), 0, pitchY);

	h_tot_allHitsCenter = new TH2D("",
			";Y-Position in Pixel [#mum];Cluster ToT", (int) (0.5 * pitchY), 0,
			pitchY, maxToT, 0, maxToT);
	h_tot_nHits1Center = new TH2D("",
			";Y-Position in Pixel [#mum];Center Pixel ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);
	h_tot_nHits2Center = new TH2D("",
			";Y-Position in Pixel [#mum];Adjacent Pixel ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);
	h_tot_nHits3PlusCenter = new TH2D("",
			";Y-Position in Pixel [#mum];Adjacent Pixels ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);

	h_sumtot_allHitsCenter = new TH2D("", ";Y-Position in Pixel [#mum];ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);
	h_sumtot_nHits1Center = new TH2D("", ";Y-Position in Pixel [#mum];ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);
	h_sumtot_nHits2Center = new TH2D("", ";Y-Position in Pixel [#mum];ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);
	h_sumtot_nHits3PlusCenter = new TH2D("", ";Y-Position in Pixel [#mum];ToT",
			(int) (0.5 * pitchY), 0, pitchY, maxToT, 0, maxToT);

	h_fractot_nHits2Center = new TH2D("",
			";Y-Position in Pixel [#mum];Fraction of ToT in Neighboring Pixel",
			(int) (0.5 * pitchY), 0, pitchY, 100, 0, 1);

}

void QShare1D::event(const TbConfig &config, const Event &event) {

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
	// Check that we have only 1 cluster:
	if (event.clusters.size() != 1) {
		return;
	}
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);
	if (matched == -1) {
		return;
	}
	// Require  hits, but apparently not too many
	if (event.hits.size() == 0 or event.hits.size() > 4) {
		return;
	}

	// Transpose so that (0,0) is bottom-left corner
	// of bottom-left edge pixel (truncated to normal size)
	double pitchX = event.dut->getPitchX();
	double pitchY = event.dut->getPitchY();
	double trackXMod = event.trackX + 0.5 * pitchX;
	double trackYMod = event.trackY + 0.5 * pitchY;
	while (trackXMod > pitchX) {
		trackXMod -= pitchX;
	}
	while (trackYMod > pitchY) {
		trackYMod -= pitchY;
	}
	double edgeDistY = trackYMod;
	double edgeDistX = trackXMod;
	// old hard-coded edge x-distance cut values
	//if( edgeDistX < 50.00 or edgeDistX > 350.00 ) { return; }
	// new cuts which maintain same percentage of x-distance removed for fe-i3 and fe-i4
	// this is a temporary hack until i know better (BJD)
	if (edgeDistX < (0.125 * pitchX) or edgeDistX > (0.875 * pitchX)) {
		return;
	}

	// Loop over all hits:
	int hitCount = 0;
	int sumtot(0), totCenter(0), totAbove(0), totBelow(0);
	bool center(false), above(false), below(false);

	for (int ii = 0; ii < event.hits.size(); ii++) {
		hitCount++;
		sumtot += event.hits.at(ii)->tot;
		double diffX = event.hits.at(ii)->col * pitchX - event.trackX;
		double diffY = event.hits.at(ii)->row * pitchY - event.trackY;
		bool xCenter(false), yCenter(false), yAbove(false), yBelow(false);

		if (fabs(diffX) < 0.5 * pitchX) {
			xCenter = true;
		}
		if (fabs(diffY) < 0.5 * pitchY) {
			yCenter = true;
		}
		if (diffY > 0.5 * pitchY and diffY < 1.5 * pitchY) {
			yAbove = true;
		}
		if (diffY < -0.5 * pitchY and diffY > -1.5 * pitchY) {
			yBelow = true;
		}

		if (xCenter and yCenter) {
			center = true;
			totCenter = event.hits.at(ii)->tot;
		}
		if (xCenter and yAbove) {
			above = true;
			totAbove = event.hits.at(ii)->tot;
		}
		if (xCenter and yBelow) {
			below = true;
			totBelow = event.hits.at(ii)->tot;
		}
	}

	h_tracks->Fill(edgeDistY);
	if (hitCount > 0 and center) {
		h_allHitsCenter->Fill(edgeDistY);
		h_sumtot_allHitsCenter->Fill(edgeDistY, sumtot);
		h_tot_allHitsCenter->Fill(edgeDistY, sumtot);
	}
	if (hitCount == 1 and center) {
		h_nHits1Center->Fill(edgeDistY);
		h_sumtot_nHits1Center->Fill(edgeDistY, sumtot);
		h_tot_nHits1Center->Fill(edgeDistY, totCenter);
	}
	if (hitCount == 2 and center) {
		if (trackYMod < 0.5 * pitchY and below) {
			h_nHits2Center->Fill(edgeDistY);
			h_sumtot_nHits2Center->Fill(edgeDistY, sumtot);
			h_tot_nHits2Center->Fill(edgeDistY, totBelow);
			h_fractot_nHits2Center->Fill(edgeDistY,
					double(totBelow) / double(sumtot));
		}
		if (trackYMod >= 0.5 * pitchY and above) {
			h_nHits2Center->Fill(edgeDistY);
			h_sumtot_nHits2Center->Fill(edgeDistY, sumtot);
			h_tot_nHits2Center->Fill(edgeDistY, totAbove);
			h_fractot_nHits2Center->Fill(edgeDistY,
					double(totAbove) / double(sumtot));
		}
	}
	if (hitCount > 2 and center) {
		if (trackYMod < 0.5 * pitchY and below) {
			h_nHits3PlusCenter->Fill(edgeDistY);
			h_sumtot_nHits3PlusCenter->Fill(edgeDistY, sumtot);
			h_tot_nHits3PlusCenter->Fill(edgeDistY, totBelow);
		}
		if (trackYMod >= 0.5 * pitchY and above) {
			h_nHits3PlusCenter->Fill(edgeDistY);
			h_sumtot_nHits3PlusCenter->Fill(edgeDistY, sumtot);
			h_tot_nHits3PlusCenter->Fill(edgeDistY, totAbove);
		}
	}
	if (hitCount > 0 and !center) {
		h_allHitsNotCenter->Fill(edgeDistY);
	}

}

void QShare1D::finalize(const TbConfig &config) {

	// Get overall charge sharing probability
	double nAllHits = h_allHitsCenter->GetEntries();
	double nSharedHits = h_nHits2Center->GetEntries()
			+ h_nHits3PlusCenter->GetEntries();

	TBALOG(kINFO) << "Tracks:                 " << h_tracks->GetEntries()
			<< endl;
	TBALOG(kINFO) << "All Hits:               " << nAllHits << endl;
	TBALOG(kINFO) << "Shared Hits:            " << nSharedHits << endl;
	TBALOG(kINFO) << "Overall Charge Sharing: " << nSharedHits / nAllHits
			<< endl;

	h_nHits2Center->Add(h_nHits3PlusCenter);
	h_nHits2Center->Divide(h_allHitsCenter);
	h_nHits2Center->SetTitle(
			";Y-Position in Pixel [#mum];Charge Sharing Fraction");
	h_nHits2Center->SetMinimum(0.0);
	h_nHits2Center->SetMaximum(1.0);
	config.drawAndSave(name, "chargeSharingCenter", h_nHits2Center);

	config.drawAndSave(name, "edgeDistY", h_tracks);

	config.drawAndSave(name, "tot_allHitsCenter", h_tot_allHitsCenter);
	config.drawAndSave(name, "tot_nHits1Center", h_tot_nHits1Center);
	config.drawAndSave(name, "tot_nHits2Center", h_tot_nHits2Center);
	config.drawAndSave(name, "tot_nHits3PlusCenter", h_tot_nHits3PlusCenter);

	config.drawAndSave(name, "sumTot_allHitsCenter", h_sumtot_allHitsCenter);
	config.drawAndSave(name, "sumTot_nHits1Center", h_sumtot_nHits1Center);
	config.drawAndSave(name, "sumTot_nHits2Center", h_sumtot_nHits2Center);
	config.drawAndSave(name, "sumTot_nHits3PlusCenter",
			h_sumtot_nHits3PlusCenter);

	config.drawAndSave(name, "fracTot_nHits2Center", "colz",
			h_fractot_nHits2Center);

	drawTotSharing(config, name, "qShare_nHits2", "", h_sumtot_nHits2Center,
			h_tot_nHits2Center, h_fractot_nHits2Center);

}

void QShare1D::drawTotSharing(const TbConfig& config, const char* analysisName,
		const char* histoName, const char* drawOpts, TH2D* h_sumTot,
		TH2D* h_otherTot, TH2D* h_fracTot) const {

	// draw profiles
	TCanvas* canvas = new TCanvas("", "", 700, 500);
	TProfile* pr_sumTot = h_sumTot->ProfileX("pr_sumTot");
	pr_sumTot->SetTitle(";Y-Position in Pixel [#mum];Mean ToT in Cluster");
	pr_sumTot->Draw();
	TString hname = TString::Format("%s_sumTot", histoName);
	char* fileName = config.buildHistName(analysisName, hname);
	canvas->SaveAs(fileName);
	delete canvas;
	delete fileName;
	delete pr_sumTot;

	canvas = new TCanvas("", "", 700, 500);
	TProfile* pr_otherTot = h_otherTot->ProfileX("pr_otherTot");
	pr_otherTot->SetTitle(
			";Y-Position in Pixel [#mum];Mean ToT in Adjacent Pixel");
	pr_otherTot->Draw();
	hname = TString::Format("%s_otherTot", histoName);
	fileName = config.buildHistName(analysisName, hname);
	canvas->SaveAs(fileName);
	delete canvas;
	delete fileName;
	delete pr_otherTot;

	canvas = new TCanvas("", "", 700, 500);
	TProfile* pr_fracTot = h_fracTot->ProfileX("pr_fracTot");
	pr_fracTot->SetStats(false);
	pr_fracTot->SetTitle(
			";Y-Position in Pixel [#mum];Mean Fraction of ToT in Adjacent Pixel");
	pr_fracTot->Draw();
	hname = TString::Format("%s_fracTot", histoName);
	fileName = config.buildHistName(analysisName, hname);
	canvas->SaveAs(fileName);
	delete canvas;
	delete fileName;
	delete pr_fracTot;

	//Now draw the distributions for each distance bin
	//TH1D* hdtot[h_sumTot->GetNbinsX()];
	assert(h_sumTot->GetNbinsX() == h_fracTot->GetNbinsX());
	for (int d = 1; d != h_sumTot->GetNbinsX() + 1; ++d) {
		canvas = new TCanvas("", "", 700, 500);
		TH1D* hs = (TH1D*) h_sumTot->ProjectionY(
				TString::Format("%i_sumtot_d=%i", this->iden, d), d, d, "e");
		hs->Rebin(4);
		//h->SetLineColor(d);
		hs->SetTitle(";ToT;Tracks");
		//hs->Draw("hist");
		hs->DrawNormalized("hist");

		TH1D* ho = (TH1D*) h_otherTot->ProjectionY(
				TString::Format("%i_othertot_d=%i", this->iden, d), d, d, "e");
		ho->Rebin(4);
		//ho->SetLineColor(d);
		ho->SetLineStyle(2);
		//ho->Draw("same hist");
		ho->DrawNormalized("same hist");

		TLegend* leg = new TLegend(0.78, 0.6, 0.98, 0.8);
		leg->SetFillColor(kWhite);
		leg->AddEntry(hs, "Sum ToT", "L");
		leg->AddEntry(ho, "Shared ToT", "L");
		leg->Draw();
		tbutils::myText(0.15, 0.85, kBlack,
				TString::Format("y = %.1f #mum",
						h_sumTot->GetXaxis()->GetBinCenter(d)));
		hname = TString::Format("%s_overlay%i", histoName, d);
		fileName = config.buildHistName(analysisName, hname);
		canvas->SaveAs(fileName);
		delete canvas;
		delete fileName;
		delete leg;

		canvas = new TCanvas("", "", 700, 500);
		TH1D* hf = (TH1D*) h_fracTot->ProjectionY(
				TString::Format("%i_frac_d=%i", this->iden, d), d, d, "e");
		hf->SetTitle(";Fraction of ToT;Tracks");
		hf->Rebin(4);
		hf->Draw();
		tbutils::myText(0.15, 0.85, kBlack,
				TString::Format("y = %.1f #mum",
						h_sumTot->GetXaxis()->GetBinCenter(d)));
		hname = TString::Format("%s_frac%i", histoName, d);
		fileName = config.buildHistName(analysisName, hname);
		canvas->SaveAs(fileName);
		delete canvas;
		delete fileName;

		delete hs;
		delete ho;
		delete hf;
	}

}
