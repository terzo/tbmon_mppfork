#ifndef ETAWIDTH_H
#define ETAWIDTH_H

// standard header files
#include <cmath>
#include <vector>

// root header files
#include <TH1D.h>
#include <TGraph.h>

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Unknown - check!
 *
 */
class EtaWidth: public TbAnalysis {

private:
	vector<TH1D*> residualsY2c;
	vector<TH1D*> residualsY2b;
	vector<TH1D*> residualsY3b;
	vector<TH1D*> residualsY3c;
	TH1D* residualY;
	TH1D* residualY22;
	TH1D* etaY;
	double centralCorrection(double ecorr, double width);
	double borderCorrection(double ecorr, double width);
	void initHisto(vector<TH1D*>& histoc, vector<TH1D*>& histob,
			double halfWidth);
	void getStats(const TbConfig& config, vector<TH1D*> h, const char* name);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //ETAWIDTH_H
