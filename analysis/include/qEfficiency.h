/*
 * qEfficiency.h
 *
 *  Created on: 16. juni 2010
 *      Author: kyrre
 */

#ifndef QEFFICIENCY_H
#define QEFFICIENCY_H

// standard header files
#include <cmath>
#include <vector>

// root header files
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "TMath.h"
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TPaletteAxis.h>

// tbmon header files
#include "clusters.h"
#include "tbanalysis.h"

/*! \brief Plots space resolved ToT and charge maps.
 *
 */
class qEfficiency: public TbAnalysis {
public:
	virtual void init(TbConfig& config);
	virtual void event(const TbConfig& config, const Event& event);
	virtual void finalize(const TbConfig& config);
private:
	TH3D* h_totScatter;
	TH3D* h_chargeScatter;
	TH1D* h_chargeLandau;

	TH3D* h_totScatterBiasLow;
	TH3D* h_totScatterBiasHi;
	TH3D* h_chargeScatterBiasLow;
	TH3D* h_chargeScatterBiasHi;

	TH3D* h_totScatterReadout;
	TH3D* h_chargeScatterReadout;

	TH2D* h_xyColumns;
	TH2D* h_xyBulk;

	//Histogram resolution
	double xRes, yRes, totRes, qRes;
	//Other config
	double electrodeWidth;
	double xy_totCutoff;

	//Common functions done oh-so-many times
	void fixTitles(TH3D* h, const char* Ztitle) {
		h->SetXTitle("Xmod position [#mu m]");
		h->SetYTitle("Ymod position [#mu m]");
		h->SetZTitle(Ztitle);
	}
	void zxProject(TH3D* h, const char* saveName, const TbConfig& config,
			TH3D* h2 = NULL);

	void drawAspect(TH2D* h, const char* saveName, const TbConfig& config,
			const Event &event);

	TProfile* makeLandauProfile(const TbConfig& config, TH3D* h_scatter,
			int rebin = 1);
};

Double_t langaufunction(Double_t *x, Double_t *par);
#endif //QEFFICIENCY_H
