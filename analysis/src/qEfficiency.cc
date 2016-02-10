/*
 * qEfficiency.cc
 *
 *  Created on: 16. June 2010, 20. Aug 2010
 *      Author: Kyrre Ness Sjøbæk
 */
#include "qEfficiency.h"

using namespace std;

void qEfficiency::init(TbConfig& config) {
	//Get -P:key val cmd line arguments
	char key[100];

	sprintf(key, "A_qEfficiency_electrodeWidth_%u", this->iden);
	electrodeWidth = config.cmdLineExtras_argGetter(key, 5.0,
			"Half-width of strip to use for profile plot [um]");

	sprintf(key, "A_qEfficiency_xy_totCutoff_%u", this->iden);
	xy_totCutoff = config.cmdLineExtras_argGetter(key, 25.0,
			"Divider for which ToT-values belong to columns or bulk");

	//Make histograms
	const DUT* dut = config.getDut(this->iden);
	this->xRes = this->yRes = 4;
	this->totRes = 1;
	this->qRes = 750;

	h_totScatter = new TH3D("h_totScatter",
			";Long Pixel Pitch [#mum];Short Pixel Pitch [#mum];TOT",
			(int) (dut->getPitchX() / this->xRes), 0, dut->getPitchX(),
			(int) (dut->getPitchY() / this->yRes), 0, dut->getPitchY(),
			(int) (150 / this->totRes), 0, 150);
	fixTitles(h_totScatter, "ToT");

	h_chargeScatter = new TH3D("h_chargeScatter",
			";Long Pixel Pitch [#mum];Short Pixel Pitch [#mum];Charge",
			(int) (dut->getPitchX() / this->xRes), 0, dut->getPitchX(),
			(int) (dut->getPitchY() / this->yRes), 0, dut->getPitchY(),
			(int) (60e3 / this->qRes), 0, 60e3);
	fixTitles(h_chargeScatter, "Charge [electrons]");

	h_chargeLandau = new TH1D("h_chargeLandau", ";Charge [ke];N",
			(int) (60e3 / this->qRes), 0, 60);
	/* FIXME
	 h_totScatterBiasLow    = new TH3D("","", (int)(400/xRes),0,400, (int)(electrodeWidth/yRes),0,electrodeWidth,                 (int)(150/totRes),0,150);
	 fixTitles(h_totScatterBiasLow,"ToT");
	 h_totScatterBiasHi     = new TH3D("","", (int)(400/xRes),0,400, (int)(electrodeWidth/yRes),0,electrodeWidth,                 (int)(150/totRes),0,150);
	 fixTitles(h_totScatterBiasHi,"ToT");
	 h_totScatterReadout    = new TH3D("","", (int)(400/xRes),0,400, (int)(2*electrodeWidth/yRes),-electrodeWidth,electrodeWidth, (int)(150/totRes),0,150);
	 fixTitles(h_totScatterReadout,"ToT");
	 h_chargeScatterBiasLow = new TH3D("","", (int)(400/xRes),0,400, (int)(electrodeWidth/yRes),50-electrodeWidth,50,             (int)(60e3/qRes),0,60e3);
	 fixTitles(h_chargeScatterBiasLow,"Charge [electrons]");
	 h_chargeScatterBiasHi  = new TH3D("","", (int)(400/xRes),0,400, (int)(electrodeWidth/yRes),50-electrodeWidth,50,             (int)(60e3/qRes),0,60e3);
	 fixTitles(h_chargeScatterBiasHi,"Charge [electrons]");
	 h_chargeScatterReadout = new TH3D("","", (int)(400/xRes),0,400, (int)(2*electrodeWidth/yRes),-electrodeWidth,electrodeWidth, (int)(60e3/qRes),0,60e3);
	 fixTitles(h_chargeScatterReadout,"Charge [electrons]");
	 h_xyColumns = new TH2D("","", (int)(400/xRes),0,400, (int)(50/yRes),0,50);
	 h_xyColumns->SetXTitle("tbutils::getFoldedX(event) position [#mu m]");
	 h_xyColumns->SetXTitle("tbutils::getPixelY(event) position [#mu m]");
	 h_xyBulk    = new TH2D("","", (int)(400/xRes),0,400, (int)(50/yRes),0,50);
	 h_xyBulk   ->SetXTitle("tbutils::getFoldedX(event) position [#mu m]");
	 h_xyBulk   ->SetXTitle("tbutils::getPixelY(event) position [#mu m]");
	 */
}

void qEfficiency::event(const TbConfig& config, const Event& event) {
	if (event.fTrack != event::kGood)
		return;
	if (event.fClusters != event::kGood)
		return;
	if (event.hits.size() == 0)
		return;
	if (event.fTrackCentralRegion != event::kGood)
		return;
	if (event.fTrackRegion != event::kGood)
		return;

	int clusterIndex = cluster::getMatched(event.clusters, event);
	if (clusterIndex == -1)
		return;

	//Fill full-pixel histos
	h_totScatter->Fill(tbutils::getFoldedX(event), tbutils::getPixelY(event),
			cluster::getSumTot(event.clusters[clusterIndex]));

	if (event.dut->getToTcalib()->hasPerPixelCalib()) {
		h_chargeScatter->Fill(tbutils::getFoldedX(event),
				tbutils::getPixelY(event),
				cluster::getSumChargePP(event.clusters[clusterIndex], event));
		h_chargeLandau->Fill(
				cluster::getSumChargePP(event.clusters[clusterIndex], event)
						/ 1000.);
	} else if (event.dut->getToTcalib()->hasCalib()) {
		h_chargeScatter->Fill(tbutils::getFoldedX(event),
				tbutils::getPixelY(event),
				cluster::getSumCharge(event.clusters[clusterIndex], event));
		h_chargeLandau->Fill(
				cluster::getSumCharge(event.clusters[clusterIndex], event)
						/ 1000.);
	}
	/* FIXME
	 //Fill partial histos
	 if (tbutils::getPixelY(event) < electrodeWidth) {
	 //Low bias row
	 h_totScatterBiasLow       ->Fill( tbutils::getFoldedX(event),tbutils::getPixelY(event), cluster::getSumTot   (event.clusters[clusterIndex]) );
	 if (event.dut->roottotcalib->hasCalib()){
	 h_chargeScatterBiasLow  ->Fill( tbutils::getFoldedX(event),tbutils::getPixelY(event), cluster::getSumCharge(event.clusters[clusterIndex],event));
	 }
	 }
	 else if ( tbutils::getPixelY(event) > (25-electrodeWidth) && tbutils::getPixelY(event) < (25+electrodeWidth) ) {
	 //Middle readout row
	 h_totScatterReadout     ->Fill( tbutils::getFoldedX(event),tbutils::getPixelY(event), cluster::getSumTot   (event.clusters[clusterIndex]) );
	 if (event.dut->roottotcalib->hasCalib()){
	 h_chargeScatterReadout->Fill( tbutils::getFoldedX(event),tbutils::getPixelY(event), cluster::getSumCharge(event.clusters[clusterIndex],event));
	 }
	 }
	 else if( tbutils::getPixelY(event) > (50-electrodeWidth) ) {
	 //Middle readout row
	 h_totScatterBiasHi      ->Fill( tbutils::getFoldedX(event),tbutils::getPixelY(event), cluster::getSumTot   (event.clusters[clusterIndex]) );
	 if (event.dut->roottotcalib->hasCalib()){
	 h_chargeScatterBiasHi ->Fill( tbutils::getFoldedX(event),tbutils::getPixelY(event), cluster::getSumCharge(event.clusters[clusterIndex],event));
	 }
	 }

	 //Fill column-histos
	 if (cluster::getSumTot(event.clusters[clusterIndex]) > xy_totCutoff) {
	 h_xyBulk->Fill(tbutils::getFoldedX(event),tbutils::getPixelY(event));
	 }
	 else {
	 h_xyColumns->Fill(tbutils::getFoldedX(event),tbutils::getPixelY(event));
	 }
	 */
}

void qEfficiency::finalize(const TbConfig& config) {
	//Save raw data to file
	config.saveToFile(this->name, "totScatter", h_totScatter);
	config.saveToFile(this->name, "chargeScatter", h_chargeScatter);
	config.saveToFile(this->name, "chargeLandau", h_chargeLandau);

	TF1 *ffit = new TF1("Langaus", langaufunction,
			(h_chargeLandau->GetMean()) / 2., 40, 4);
	Double_t startvalues[4] =
			{ 0.65, h_chargeLandau->GetMean(), h_chargeLandau->Integral(
					(h_chargeLandau->GetMean()) / 2., 40), 1.6 };
	ffit->SetParameters(startvalues);
	ffit->SetParNames("Width", "MPV", "Area", "GSigma");
	h_chargeLandau->Fit("Langaus", "LLR");
	config.drawAndSave(this->name, "chargeLandau_fitted", h_chargeLandau);
	TBALOG(kINFO) << "Langaus MPV:    " << ffit->GetParameter(1) << " +- "
			<< ffit->GetParError(1) << " ke\n";
	TBALOG(kINFO) << "Langaus Width:  " << ffit->GetParameter(0) << " +- "
			<< ffit->GetParError(0) << "\n";
	TBALOG(kINFO) << "Langaus Area:   " << ffit->GetParameter(2) << " +- "
			<< ffit->GetParError(2) << "\n";
	TBALOG(kINFO) << "Langaus GSigna: " << ffit->GetParameter(3) << " +- "
			<< ffit->GetParError(3) << "\n";
	/* FIXME
	 config.saveToFile(this->name, "totScatterBiasLow"      , h_totScatterBiasLow);
	 config.saveToFile(this->name, "totScatterBiasHi"       , h_totScatterBiasHi);
	 config.saveToFile(this->name, "totScatterReadout"      , h_totScatterReadout);
	 config.saveToFile(this->name, "chargeScatterBiasHi"    , h_chargeScatterBiasHi);
	 config.saveToFile(this->name, "chargeScatterReadout"   , h_chargeScatterReadout);
	 config.saveToFile(this->name, "xyBulk"                 , h_xyBulk);
	 config.saveToFile(this->name, "xyColumns"              , h_xyColumns);
	 */

	//Project out y-axis to make x-z scatterplot and
	//make normal (mean/spread) 2D and 1D profiles
	zxProject(h_totScatter, "tot", config);
	zxProject(h_chargeScatter, "charge", config);
	/* FIXME
	 zxProject(h_totScatterBiasLow    , "totBiasLow"         , config);
	 zxProject(h_totScatterBiasHi     , "totBiasHi"          , config);
	 zxProject(h_totScatterReadout    , "totReadout"         , config);
	 zxProject(h_chargeScatterBiasLow , "chargeBiasLow"      , config);
	 zxProject(h_chargeScatterBiasHi  , "chargeBiasHi"       , config);
	 zxProject(h_chargeScatterReadout , "chargeReadout"      , config);
	 zxProject(h_totScatterBiasLow    , "totBiasCombined"    , config, h_totScatterBiasHi);
	 zxProject(h_chargeScatterBiasLow , "chargeBiasCombined" , config, h_chargeScatterBiasHi);
	 */

	//Column finders
	/*
	 config.drawToFile(this->name,"xyBulk"                 , h_xyBulk);
	 config.drawToFile(this->name,"xyBulk_colz", "colz"    , h_xyBulk);
	 config.drawToFile(this->name,"xyColumns"              , h_xyColumns);
	 config.drawToFile(this->name,"xyColumns_colz", "colz" , h_xyColumns);
	 */
//   drawAspect(h_xyBulk,"xyBulk",config,event);
//   drawAspect(h_xyColumns,"xyColumns",config,event);

	//Project out y-axis to make x-z scatterplot
	TH2D* h_totProject = (TH2D*) h_totScatter->Project3D("zx");
	TH2D* h_chargeProject = (TH2D*) h_chargeScatter->Project3D("zx");
	h_totProject->SetTitle("");
	h_chargeProject->SetTitle("");
	config.drawAndSave(this->name, "totProject", h_totProject);
	config.drawAndSave(this->name, "chargeProject", h_chargeProject);

	//Make normal (mean/spread) 2D and 1D profiles
	TProfile2D* p_totProfile = h_totScatter->Project3DProfile("yx");
	config.saveToFile(this->name, "totProfile", p_totProfile);
	config.drawToFile(this->name, "totProfile", "colz", p_totProfile);

	// BJD: this is a dirty hack for IBL sensor review request, 4 July 2011
	// Matthias George made me do it!!!
	TH2D* h2_totProfile = (TH2D*) p_totProfile->ProjectionXY();
	h2_totProfile->SetTitle(
			";Long Pixel Pitch [#mum];Short Pixel Pitch [#mum];TOT");
	const DUT* dut = config.getDut(this->iden);
	double k = dut->getPitchY() / dut->getPitchX();
	TCanvas* c_aspect = new TCanvas("", "", 1200, (int) floor(1200 * k));
	c_aspect->cd();
	gStyle->SetPalette(1, 0);
	h2_totProfile->SetStats(false);
	h2_totProfile->GetYaxis()->SetTickLength(0.01);
	// Tinker with presentation, since ROOT tries so hard to make these plots unreadable
	gPad->SetRightMargin(0.7 * k);
	gPad->SetLeftMargin(0.5 * k);
	gPad->SetBottomMargin(0.22);
	h2_totProfile->GetXaxis()->SetNdivisions((int) dut->getPitchX() / 50.0);
	h2_totProfile->GetYaxis()->SetNdivisions(5);
	h2_totProfile->GetXaxis()->SetLabelOffset(0.005);
	h2_totProfile->GetYaxis()->SetLabelOffset(0.005);
	h2_totProfile->GetZaxis()->SetLabelOffset(0.005);
	h2_totProfile->GetXaxis()->SetLabelSize(0.1);
	h2_totProfile->GetYaxis()->SetLabelSize(0.1);
	h2_totProfile->GetZaxis()->SetLabelSize(0.08);
	h2_totProfile->GetXaxis()->SetTitleOffset(1.0);
	h2_totProfile->GetYaxis()->SetTitleOffset(2.0 * k);
	h2_totProfile->GetZaxis()->SetTitleOffset(2.0 * k);
	h2_totProfile->GetXaxis()->SetTitleSize(0.1);
	h2_totProfile->GetYaxis()->SetTitleSize(0.07);
	h2_totProfile->GetZaxis()->SetTitleSize(0.07);
	h2_totProfile->Draw("colz");
	gPad->Update();
	char* fileName = config.buildHistName(this->name, "totProfileAspect");
	c_aspect->SaveAs(fileName);
	delete c_aspect;

	TProfile2D* p_chargeProfile = h_chargeScatter->Project3DProfile("yx");
	config.saveToFile(this->name, "chargeProfile", p_chargeProfile);
	config.drawToFile(this->name, "chargeProfile", "colz", p_chargeProfile);

	TProfile* p_totProfile_all = h_totProject->ProfileX();
	TProfile* p_chargeProfile_all = h_chargeProject->ProfileX();
	config.drawAndSave(this->name, "totProfile_all", p_totProfile_all);
	config.drawAndSave(this->name, "chargeProfile_all", p_chargeProfile_all);

	/*
	 //Sliced raw 1D profiles (get slices from 2D)
	 h_chargeScatter->GetYaxis()->SetRangeUser(25-this->electrodeWidth, 25+this->electrodeWidth);
	 TH2D* h_chargeSliceReadout = (TH2D*) h_chargeScatter->Project3D("zx");
	 TProfile* p_chargeSliceReadout_profile = h_chargeSliceReadout->ProfileX("",1,-1,"s");
	 config.drawAndSave(this->name,"chargeSliceReadout_profile",p_chargeSliceReadout_profile);

	 h_chargeScatter->GetYaxis()->SetRangeUser(0,this->electrodeWidth);
	 TH2D* h_chargeSliceBias1 = (TH2D*) h_chargeScatter->Project3D("zx");
	 h_chargeScatter->GetYaxis()->SetRangeUser(50-this->electrodeWidth,50);
	 TH2D* h_chargeSliceBias2 = (TH2D*) h_chargeScatter->Project3D("zx");
	 TH2D* h_chargeSliceBias = (TH2D*) h_chargeSliceBias1->Clone();
	 h_chargeSliceBias->Add(h_chargeSliceBias2);
	 TProfile* p_chargeSliceBias_profile = h_chargeSliceBias->ProfileX("",1,-1,"s");
	 config.drawAndSave(this->name,"chargeSliceBias_profile",p_chargeSliceBias_profile);
	 //Reset axis
	 h_chargeScatter->GetYaxis()->SetRange(1,h_chargeScatter->GetYaxis()->GetNbins());

	 //Make "landaufit" 2D and 1D profiles: Slice the histo into x-axis ranges
	 if (!config.isSimulation) {
	 config.drawAndSave(this->name, "chargeLandaufitProfile", makeLandauProfile(config, h_chargeScatter));
	 }

	 //Sliced landaufit 1D
	 */
}

void qEfficiency::zxProject(TH3D* h, const char* saveName,
		const TbConfig& config, TH3D* h2) {
	char buffer[100];

	//<Norwegian joke>
	//Det skakk'e værra lett, for pokker! Da kunnjo kæmsomhælst læri sæ ROOT og blitt partikkefysiker!
	//Slikt skavikke hanokkå av! Fysjom! Tænk på dom sorte hulla!
	//</Norwegian joke>
	double range_ylow = h->GetZaxis()->GetBinLowEdge(h->GetZaxis()->GetFirst());
	double range_yhi = h->GetZaxis()->GetBinUpEdge(h->GetZaxis()->GetLast());

	TH2D* h_project = (TH2D*) h->Project3D("zx");
	h_project->SetTitle("");
	if (h2 != NULL) {
		TH2D* h2_project = (TH2D*) h2->Project3D("zx");
		h_project->Add(h2_project);
	}
	sprintf(buffer, "%s_%s", saveName, "project");
	h_project->GetYaxis()->SetRangeUser(range_ylow, range_yhi);
	config.drawAndSave(this->name, buffer, h_project);

	TProfile* p_profile = h_project->ProfileX();
	sprintf(buffer, "%s_%s", saveName, "profile");
	p_profile->GetYaxis()->SetRangeUser(range_ylow, range_yhi);
	config.drawAndSave(this->name, buffer, p_profile);
}

void qEfficiency::drawAspect(TH2D* h, const char* saveName,
		const TbConfig& config, const Event &event) {
	double k = event.dut->getPitchY() / event.dut->getPitchX();
	TCanvas* c_aspect = new TCanvas("", "", 1600, (int) floor(1600 * k));
	c_aspect->cd();

	gStyle->SetPalette(1, 0);
	h->SetStats(false);

	//Remove the ticks on the axis
	//h->GetXaxis()->SetTickLength(0.00000000001);
	h->GetYaxis()->SetTickLength(0.01);
	//Reduce the axis label entries
	h->GetXaxis()->SetNdivisions(8);
	h->GetYaxis()->SetNdivisions(3);
	h->GetXaxis()->SetLabelOffset(0.005);
	h->GetYaxis()->SetLabelOffset(0.009);
	h->GetXaxis()->SetLabelSize(0.12);
	h->GetYaxis()->SetLabelSize(0.16);
	//Remove the titles
	h->SetTitle(";;");

	h->Draw("cont4z");
	gPad->Update();

	TPaletteAxis *palette = (TPaletteAxis*) h->GetListOfFunctions()->FindObject(
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

	char* fileName = config.buildHistName(this->name, saveName);
	c_aspect->SaveAs(fileName);
	delete c_aspect;

	//This is a canvas with wrong aspect ratio used to get the large palette
	TCanvas* canvas2 = new TCanvas("", "", 1600, 1600);
	canvas2->cd();
	gStyle->SetPalette(1, 0);
	//h->Draw("cont4z");
	//gPad->Update();
	palette = (TPaletteAxis*) h->GetListOfFunctions()->FindObject("palette");
	if (palette != 0) {
		palette->SetY2NDC(0.99);
		palette->SetY1NDC(0.03);
		palette->SetX2NDC(0.7);
		palette->Draw();
	} else {
		TBALOG( kERROR ) << "palette not found!" << endl;
		//exit(-1);
	}
	char* histoName2 = new char[600];
	sprintf(histoName2, "%s_palette", saveName);
	fileName = config.buildHistName(this->name, histoName2);
	canvas2->SaveAs(fileName);
	delete fileName;
	delete canvas2;

}

TProfile* qEfficiency::makeLandauProfile(const TbConfig& config,
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

	return p_landauProfile;
}

Double_t langaufunction(Double_t *x, Double_t *par) {

	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation), 
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.

	// Numeric constants
	Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
	Double_t mpshift = -0.22278298; // Landau maximum location

	// Control constants
	Double_t np = 100.0; // number of convolution steps
	Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas

	// Variables
	Double_t xx;
	Double_t mpc;
	Double_t fland;
	Double_t sum = 0.0;
	Double_t xlow, xupp;
	Double_t step;
	Double_t i;

	// MP shift correction
	mpc = par[1] - mpshift * par[0];

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp - xlow) / np;

	// Convolution integral of Landau and Gaussian by sum
	for (i = 1.0; i <= np / 2; i++) {
		xx = xlow + (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);

		xx = xupp - (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);
	}

	return (par[2] * step * sum * invsq2pi / par[3]);
}

