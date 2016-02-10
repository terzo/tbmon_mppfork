#ifndef EFFICIENCY_H
#define EFFICIENCY_H

// standard header files
#include <assert.h>

// root header files
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include <TH2D.h>
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Calculates the matched and unmatched efficiencies and other track based parameters. Plots the corresponding plots.
 *
 *
 */
class Efficiency: public TbAnalysis {

private:
	int ntracks; ///< total number of tracks in all runs
	int nhits; ///< total number of matched tracks in all runs
	int anyhit; ///< total number of all (whether matched or unmatched) hits in all runs

	TH2D* h_trackMap; ///<
	TH2D* h_hitMap; ///<
	TH2D* h_anyHitMap; ///<
	TH2D* h_trackMap216; ///<
	TH2D* h_hitMap216; ///<
	TH2D* h_anyHitMap216; ///<
	TH2D* h_trackPixelMap; ///<
	TH2D* h_hitPixelMap; ///<
	TH2D* h_trackPixelBiasMap; ///<
	TH2D* h_hitPixelBiasMap; ///<
	TH2D* h_trackPixelReadoutMap; ///<
	TH2D* h_hitPixelReadoutMap; ///<

	double effPlotMinY; ///<
	int pixMapBinSizeX; ///<
	int pixMapBinSizeY; ///<

	//
	void drawToFileSensorAspect(const TbConfig& config,
			const char* analysisName, const char* histoName,
			const char* drawOpts, TH1* h1, TH1* h2 = NULL, TH1* h3 = NULL,
			TH1* h4 = NULL) const;
	//
	void drawToFilePixelMapEff(const TbConfig& config, const char* analysisName,
			const char* histoName, TH2D* hitPixelMap, TH2D* trackPixelMap);
	//
	bool biasElecRow(const double& ypos, const Event &event);
	//
	bool readoutElecRow(const double& ypos, const Event &event);
	//
	void drawToFileElecEffComp(const TbConfig& config, const char* analysisName,
			const char* histoName, TH2D* hitPixelReadMap,
			TH2D* trackPixelReadMap, TH2D* hitPixelBiasMap,
			TH2D* trackPixelBiasMap);
	//
	void SetAspectStyle(TH1* h);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //EFFICIENCY_H
