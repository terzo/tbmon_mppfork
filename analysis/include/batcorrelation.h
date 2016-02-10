#ifndef BATCORRELATION_H
#define BATCORRELATION_H

// standard header files
#include <cmath>
#include <fstream>
#include <iostream>

// root header files
#include "TH2D.h"
#include <TProfile.h>

// tbmon header files
#include "event.h"
#include "tbanalysis.h"

/*! \brief ONLY FOR BAT DATA (no longer maintained - will be probably removed at some point): correlation plots
 *
 */
class BatCorrelation: public TbAnalysis {
private:
	//Data storage
	TProfile* h_13x;
	TProfile* h_16x;
	TProfile* h_36x;
	TProfile* h_13y;
	TProfile* h_16y;
	TProfile* h_36y;
	//Plots for plotting
	TH2D* h_16x_scatter;
	TH2D* h_16y_scatter;
	ofstream stream;
public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
	virtual void initRun(const TbConfig &config);
	virtual void finalizeRun(const TbConfig &config);
};

#endif //BATCORRELATION_H
