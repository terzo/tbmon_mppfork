#ifndef EDGEEFFICIENCYSHIFT_H
#define EDGEEFFICIENCYSHIFT_H

// standard header files
#include <assert.h>

// root header files
#include "TCanvas.h"
#include "TF1.h"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "TGraphAsymmErrors.h"
#include <TH2D.h>
#include <TH3D.h>
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include <TProfile.h>
#include <TProfile2D.h>
#include "TROOT.h"
#include "TStyle.h"

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Specialized efficiency class which plots global and sub-pixel resolved efficiency plots for
 * "pixel shifted stepwise" FE-I3 sized sensors. Those sensors feature groups of 10 pixels where the guard rings
 * are shifted below those pixels to a differing degree. For sub-pixel resolved plots 8 of those pixels
 * per group are overlayed. 2 are omitted in order to avoid peculiar field bending on the borders of
 * the groups. Does not make sense for any other kind of sensor.
 *
 */
class EdgeEfficiencyShift: public TbAnalysis {
private:
	int ntracks;
	int ntracksedge;
	int nhitsedge;

	TH2D* leftTrackPixelMapEdge;
	TH2D* leftHitPixelMapEdge;
	TH2D* rightTrackPixelMapEdge;
	TH2D* rightHitPixelMapEdge;
	TH2D* rightTrackOverlay;
	TH2D* rightHitOverlay;
	TH2D* hitOverlay;
	TH2D* trackOverlay;

	void drawToFileSensorAspect(const TbConfig& config,
			const char* analysisName, const char* histoName,
			const char* drawOpts, TH1* h1, TH1* h2 = NULL, TH1* h3 = NULL,
			TH1* h4 = NULL) const;
	void drawToFileEdgeEff(const TbConfig& config, const Event &event,
			const char* analysisName, const char* histoName,
			TObject* hist) const;
	void drawToFilePixelMapEff(const TbConfig& config, const Event &event,
			const char* analysisName, const char* histoName, TH2D* hitPixelMap,
			TH2D* trackPixelMap) const;

	int getEdgeMatched(const vector<event::cluster_t> &clusters,
			const Event &event);
	bool isEdgeMatched(const event::cluster_t &clu, const Event& event);

	TProfile* makeLandauProfile(const TbConfig& config, TH3D* h_scatter,
			int rebin = 1);

	//charge stuff
	TH3D* h_chargeScatterRight;
	TH3D* h_chargeScatterFolded;

	TProfile2D* cleanUp(const TH3D& chargeScatter);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //EDGEEFFICIENCYSHIFT_H
