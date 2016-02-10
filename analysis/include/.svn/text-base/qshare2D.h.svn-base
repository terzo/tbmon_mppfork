#ifndef QSHARE2D_H
#define QSHARE2D_H

// standard header files
#include <iostream>
#include <string>

// root header files
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TProfile2D.h"
#include "TStyle.h"

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief QShare2D produces 3 2d histograms to represent charge sharing
 *
 * All histograms are plotted with respect to track coordinates, modulo 2 by 2 pixel size
 * 
 * h_anyHitTrkMap: entry for every time the (from track coordinates) corresponding OR any adjacent pixel fires
 * h_sharedHitTrkMap: entry for every time the corresponding pixels fires AND any other adjacent pixel fires (i.e. when charge is shared)
 * The third histogram is the ratio of entries in h_sharedHitTrkMap to h_anyHitTrkMap
 *
 * The analysis regards only the events with:
 * exactly one cluster
 * between 1 and 4 hits (inclusive) in this one cluster
 * a hit in the central region away from masked pixels
 */
class QShare2D: public TbAnalysis {

private:

	TH2D* h_anyHitTrkMap;
	TH2D* h_sharedHitTrkMap;
	TH2D* h_anyHitTrkMap_mirror; //mirrored version to be aligned with eff and qeff plots
	TH2D* h_sharedHitTrkMap_mirror; //mirrored version to be aligned with eff and qeff plots
	TH3D* h_directCharge;
	TH3D* h_chargeShareScatter;

	bool doCharge;
	void drawAspect(const TbConfig& config, const char* analysisName,
			const char* histoName, TH2D* h_any, TH2D* h_shared) const;

public:

	int numTracks;

	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //QSHARE2D_H
