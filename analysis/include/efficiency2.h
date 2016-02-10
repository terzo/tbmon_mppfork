#ifndef EFFICIENCY2_H
#define EFFICIENCY2_H

#include <assert.h>

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

/*! \brief Similar to efficiency - needs to be merged!
 *
 */
class Efficiency2: public TbAnalysis {

private:
	const DUT* dut;
	double pitchX;
	double pitchY;
	int nCols;
	int nRows;

	int ntracks;
	int nhits;
	int anyhit;

	TH2D* trackMap;
	TH2D* hitMap;
	TH2D* anyhitMap;
	TH2D* trackPixelMap;
	TH2D* hitPixelMap;
	TH2D* trackPixelBiasMap;
	TH2D* hitPixelBiasMap;
	TH2D* trackPixelReadoutMap;
	TH2D* hitPixelReadoutMap;

	void drawToFileSensorAspect(const TbConfig& config, const Event &event,
			const char* analysisName, const char* histoName,
			const char* drawOpts, TH1* h1, TH1* h2 = NULL, TH1* h3 = NULL,
			TH1* h4 = NULL) const;
	void drawToFilePixelMapEff(const TbConfig& config, const char* analysisName,
			const char* histoName, TH2D* hitPixelMap, TH2D* trackPixelMap);
	bool biasElecRow(const double& ypos, const Event &event);
	bool readoutElecRow(const double& ypos);
	void drawToFileElecEffComp(const TbConfig& config, const char* analysisName,
			const char* histoName, TH2D* hitPixelReadMap,
			TH2D* trackPixelReadMap, TH2D* hitPixelBiasMap,
			TH2D* trackPixelBiasMap);
	void SetAspectStyle(TH1* h);

	double effplots_lowRange;
	int effplots_binSizeX;
	int effplots_binSizeY;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //EFFICIENCY2_H
