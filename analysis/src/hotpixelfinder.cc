#include "hotpixelfinder.h"

void HotPixelFinder::init(TbConfig &config) {

	TBALOG(kDEBUG) << "Entered hotpixel init." << std::endl;

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();
	nCols = dut->getNcols();
	nRows = dut->getNrows();

	eventCount = 0;
	numberOfHotPixels = 0;
	numberOfDeadPixels = 0;

	// get options or set standard values
	_maximumOccupancy =
			(double) config.cmdLineExtras_argGetter(
					(string) "A_hotpixelfinder_maximumOccupancy",
					(double) 0.5e-3,
					(string) "Minimum occupancy which leads to flagging a pixel as hot.");

	//Extract lv1 configuration from config
	lv1Min = dut->lv1Min;
	lv1Max = dut->lv1Max;
	TBALOG(kINFO) << "Cuts from DUT calib: " << lv1Min << " < lv1 < " << lv1Max
			<< endl;
	//Check if DUT is properly configured
	if (lv1Max < lv1Min or lv1Max == -1) {
		TBALOG(kERROR) << "Check the lv1 calibration for DUT " << iden << endl
				<< "No way this is going to work. Try again." << endl;
		exit(-1);
	}

	char title[300];

	h_hitmap = new TH2D("", ";Column;Row", nCols, 0, nCols, nRows, 0, nRows);
	h_hitmap->SetZTitle("Hits");

	h_masked = new TH2D("", ";Column;Row", nCols, 0, nCols, nRows, 0, nRows);
	sprintf(title,"0 = unmasked / 1 = dead / 2 = noisy (cut at >= %g)",_maximumOccupancy);
	h_masked->SetZTitle(title);

	h_outOfTime = new TH2D("", ";Column;Row", nCols, 0, nCols, nRows, 0, nRows);
	h_outOfTime->SetZTitle("Hits");

	h_noiseOccupancy = new TH1D("", ";;Pixels", 10000, 0, 1);
	sprintf(title,"Noise occupancy (using out of time criterium, cut at >= %g)",_maximumOccupancy);
	h_noiseOccupancy->GetXaxis()->SetTitle(title);

	char exec_lv1_text[1000];
	sprintf(exec_lv1_text,"double current_XMin = gPad->GetUxmin();double current_YMin = gPad->GetUymin();double current_XMax = gPad->GetUxmax();double current_YMax = gPad->GetUymax();double X1;double Y1;double X2;double Y2;double X1_nom = (double)%i;double X2_nom=(double)%i;if(X1_nom>current_XMin){X1=X1_nom;} else {X1=current_XMin;};if(X2_nom<current_XMax){X2=X2_nom;}else{X2=current_XMax;};Y1=current_YMin;Y2=current_YMax;if(X2>X1){if(gPad->FindObject(\"TBox\")){TBox* _inTime = (TBox*)gPad->FindObject(\"TBox\");_inTime->SetX1(X1);_inTime->SetX2(X2);_inTime->SetY1(Y1);_inTime->SetY2(Y2);} else{TBox* _inTime = new TBox(X1,Y1,X2,Y2);_inTime->SetFillStyle(3351);_inTime->SetFillColor(kGreen);_inTime->Draw(\"SAME\");}}",lv1Min,lv1Max+1); // added +1 by Botho
	//sprintf(exec_lv1_text,"double current_XMin = gPad->GetUxmin();double current_YMin = gPad->GetUymin();double current_XMax = gPad->GetUxmax();double current_YMax = gPad->GetUymax();double X1;double Y1;double X2;double Y2;double X1_nom = (double)%i;double X2_nom=(double)%i;if(X1_nom>current_XMin){X1=X1_nom;} else {X1=current_XMin;};if(X2_nom<current_XMax){X2=X2_nom;}else{X2=current_XMax;};Y1=current_YMin;Y2=current_YMax;printf(\"XY(%%f,%%f,%%f,%%f)\\n\",X1,Y1,X2,Y2);printf(\"NOMX(%%f,%%f)\\n\",X1_nom,X2_nom);printf(\"CURRENT(%%f,%%f,%%f,%%f)\\n\",current_XMin,current_YMin,current_XMax,current_YMax);if(X2>X1){if(gPad->FindObject(\"TBox\")){TBox* _inTime = (TBox*)gPad->FindObject(\"TBox\");_inTime->SetX1(X1);_inTime->SetX2(X2);_inTime->SetY1(Y1);_inTime->SetY2(Y2);} else{TBox* _inTime = new TBox(X1,Y1,X2,Y2);_inTime->SetFillStyle(3351);_inTime->SetFillColor(kGreen);_inTime->Draw(\"SAME\");}}",lv1Min,lv1Max);
	TExec* lv1_ex = new TExec("exec_lv1",exec_lv1_text);
	h_lv1 = new TH1D("", ";lv1;Hits", 16, 0, 16);
	h_lv1->GetListOfFunctions()->Add(lv1_ex);
}


void HotPixelFinder::event(const TbConfig &config, const Event &event) {

	const DUT* dut = config.getDut(this->iden);
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	eventCount++;
	for (int ii = 0; ii < event.rawHits.size(); ii++) {
		h_lv1->Fill(event.rawHits.at(ii)->lv1);
		h_hitmap->Fill(event.rawHits.at(ii)->col, event.rawHits.at(ii)->row);

		// lv1 cut changed by Botho
		// to lv1Min <= lv1 <= lv1Max
		// so that lv1Min = 0, lv1Max = 15 does not cut any events
		/*if (event.rawHits.at(ii)->lv1 > lv1Min
				and event.rawHits.at(ii)->lv1 < lv1Max) {
		} */
		if (event.rawHits.at(ii)->lv1 >= lv1Min
				and event.rawHits.at(ii)->lv1 <= lv1Max) {
		} else
			h_outOfTime->Fill(event.rawHits.at(ii)->col,
					event.rawHits.at(ii)->row);

	}

}

void HotPixelFinder::checkDead() {

	for (int xx = 1; xx <= h_hitmap->GetNbinsX(); xx++) {
		for (int yy = 1; yy <= h_hitmap->GetNbinsY(); yy++) {
			if (h_hitmap->GetBinContent(xx, yy) == 0) {
				h_masked->SetBinContent(xx, yy, 1);
				stream << xx - 1 << ":" << yy - 1 << endl;
				numberOfDeadPixels++;
			}
		}
	}

}

void HotPixelFinder::checkNoise() {

	h_outOfTime->Scale(1.0 / eventCount);
	
	// why is this scaling applied? it also introduces potential division by 0...
	// commented out by Botho
//	h_outOfTime->Scale(16.0 / (15 - (lv1Max - lv1Min)));

	for (int xx = 1; xx <= h_outOfTime->GetNbinsX(); xx++) {
		for (int yy = 1; yy <= h_outOfTime->GetNbinsY(); yy++) {
			double val = h_outOfTime->GetBinContent(xx, yy);
			if (val > 0) {
				h_noiseOccupancy->Fill(val);
			}
			if (val > _maximumOccupancy) {
				h_masked->SetBinContent(xx, yy, 2);
				stream << xx - 1 << ":" << yy - 1 << endl;
				numberOfHotPixels++;
			}
		}
	}

}

void HotPixelFinder::finalize(const TbConfig &config) {

	stream.open(config.getOutStreamName(this->name, "masks"));
	checkDead();
	checkNoise();
	stream.close();

	char* histoTitle = new char[500];
	sprintf(histoTitle, "Raw hitmap %i", this->iden);
	h_hitmap->SetTitle(histoTitle);
	config.drawAndSave(this->name, "hitmap", "colz", h_hitmap);

	sprintf(histoTitle, "Noise occupancy %i", this->iden);
	h_noiseOccupancy->SetTitle(histoTitle);
	config.drawAndSave(this->name, "1D-occupancy", h_noiseOccupancy);

	h_outOfTime->SetStats(false);
	sprintf(histoTitle, "Noise occupancy %i", this->iden);
	h_outOfTime->SetTitle(histoTitle);
	//h_outOfTime->GetZaxis()->SetRangeUser(0.0,1e-4);
	config.drawAndSave(this->name, "outoftime", "colz", h_outOfTime);

	sprintf(histoTitle, "Masks %i", this->iden);
	h_masked->SetTitle(histoTitle);
	h_masked->SetStats(false);
	h_masked->GetZaxis()->SetRangeUser(0.0, 2.0);
	config.drawAndSave(this->name, "occupancy-masks", "colz", h_masked);

	sprintf(histoTitle, "Lv1 %i", this->iden);
	h_lv1->SetTitle(histoTitle);
	config.drawAndSave(this->name, "lv1", h_lv1);

	double percentageDead = (((double) numberOfDeadPixels
			/ (double) (nCols * nRows)) * 100);
	double percentageHot =
			(((double) numberOfHotPixels / (nCols * nRows)) * 100);
	TBALOG(kINFO) << std::endl;
	TBALOG(kINFO) << "Number of pixels identified as hot: "
			<< numberOfHotPixels << "(" << percentageHot << "%)" << std::endl;
	TBALOG(kINFO) << "Number of pixels identified as dead: "
			<< numberOfDeadPixels << "(" << percentageDead << "%)" << std::endl;

}
