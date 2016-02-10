#ifndef BATUNBIASED_H
#define BATUNBIASED_H
/*
 * This analysis plots unbiased residuals for BAT,
 * useful for knowing how much faith to put into
 * simulation resolution estimates.
 *
 * Kyrre N. Sjobak (k.n.sjobak@fys.uio.no)
 */

// standard header files
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

// root header files
#include <TH1D.h>
#include <TH2D.h>

// tbmon header files
#include <tbanalysis.h>
#include "Track.h"

/*! \brief ONLY FOR BAT DATA (no longer maintained - will be probably removed at some point): plots unbiased residuals for BAT,
 * useful for knowing how much faith to put into simulation resolution estimates.
 *
 */
class BatUnbiased: public TbAnalysis {
private:
	static const int batIdens[3];
	static const int numBats = 3;

	TH1D*** h_residuals; //First dimension: plane
						 //Second dimension: layer (y,x)
	TH1D*** h_residuals_zoom;

	TH2D*** h2_residualsVsEtaCorr; //Residuals as a function
								   // of eta-corrected hit position,
								   // folded into one strip pitch
	TH1D*** h_etaCorr; //Check that eta-corrected hit position is ~flat

	//Histos that are used with simulation data
	TH1D*** h_trackVsTruth;
	TH1D*** h_hitVsTruth;

	//Plot the true residual as a function of where the particle hit ("true-eta")
	TH2D*** h2_trueresidualVSinterstrip_truehit;
	TH2D*** h2_trueresidualVsEta;

	//Plot the true residual as a function of the eta-corrected hit position
	TH2D*** h2_trueresidualVsEtaCorr;

	//Find the correlations:
	//1st: Plane
	//2nd: layer
	//3rd: x,y of truth
	TH2D**** h2_stripVstruth;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //BATUNBIASED_H
