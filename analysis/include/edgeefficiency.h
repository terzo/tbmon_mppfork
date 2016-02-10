#ifndef EDGEEFFICIENCY_H
#define EDGEEFFICIENCY_H

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

/*! \brief Global and sub-pixel resolved efficiency plots for the edge region.
 *
 */
class EdgeEfficiency: public TbAnalysis {

private:
	int nTracks;
	int nTracksEdge;
	int nHitsEdge;

	TH2D* h_EEAnyHitSensorMap;
	TH2D* h_EEHitSensorMap;
	TH2D* h_EETrackSensorMap;
	TH2D* h_EESensorMap;
	TH2D* h_EEHitPixelMap;
	TH2D* h_EEHitLeftPixelMap;
	TH2D* h_EEHitRightPixelMap;
	TH2D* h_EETrackPixelMap;
	TH2D* h_EETrackLeftPixelMap;
	TH2D* h_EETrackRightPixelMap;

	void drawToFileEdgeEff(const TbConfig& config, const char* analysisName,
			const char* histoName, TObject* hist, double fitLeft,
			double fitRight) const;

	void ApplySensorMapStyle(const TbConfig& config, TH2D* histo,
			const char* title, const char* fileName,
			const char* histoName) const;
	void ApplyPixelMapStyle(const TbConfig& config, TH2D* histo,
			const char* title, const char* fileName,
			const char* histoName) const;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //EDGEEFFICIENCY_H
