#ifndef EFFICIENCYVSRUN_H
#define EFFICIENCYVSRUN_H

// standard header files
#include <assert.h>
#include <map>

// root header files
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include <TGraphErrors.h>
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

using namespace std;

/*! \brief Plots matched and unmatched efficiency vs. run number
 *
 */
class EfficiencyVsRun: public TbAnalysis {

private:
	TH1D* histo_EffAllVsRun;
	TH1D* histo_EffMatchedVsRun;

	// n(tracks) vs run
	map<int, int> ntracks;
	// n(all hits) vs run
	map<int, int> nhits;
	// n(matched hits) vs run
	map<int, int> nmatched;

	// efficiency (all hits) vs run
	TGraphErrors* h_EffAllVsRun;
	// efficiency (matched hits) vs run
	TGraphErrors* h_EffMatchedVsRun;

	void drawAndSave(const TbConfig& config, const char* analysisName,
			const char* histoName, TGraphErrors* gr);

	double effPlotMinY;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //EFFICIENCYVSRUN_H
