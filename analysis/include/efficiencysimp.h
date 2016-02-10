#ifndef EFFICIENCYSIMP_H
#define EFFICIENCYSIMP_H

#include <assert.h>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include <TH2D.h>
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"

#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Plots charge and hit maps - redundant?
 *
 */
class EfficiencySimp: public TbAnalysis {
private:
	int ntracks;
	int nhits;
	int anyhit;

	vector<int>* totStore[18][160];
	vector<float>* chargeClusterStore[18][160];
	vector<int> lvl1;

	TH1I* totHist;
	TTree* chargeClusterTree;
	TH1I* lvl1Hist;
	TH2D* totMap;
	TH2D* chargeClusterMap;
	TH2D* trackMap;
	TH2D* hitMap;
	TH2D* anyhitMap;
	TH2D* trackPixelMap;
	TH2D* hitPixelMap;
	TH2D* trackPixelBiasMap;
	TH2D* hitPixelBiasMap;
	TH2D* trackPixelReadoutMap;
	TH2D* hitPixelReadoutMap;

	void readMask(string maskpath, TbConfig & config);

	double effplots_lowRange;
	int effplots_binSizeX;
	int effplots_binSizeY;

	bool lvl1cut;
	bool maskon;
	bool chargeconversion;
	bool **mask;
	fstream maskfile;
	fstream calibfile;

	float calA[18][160];
	float calB[18][160];
	float calC[18][160];

	int chargesum;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //EFFICIENCYSIMP_H
