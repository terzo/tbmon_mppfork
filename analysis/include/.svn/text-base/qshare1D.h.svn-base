#ifndef QSHARE1D_H
#define QSHARE1D_H

// standard header files
#include <assert.h>
#include <iostream>
#include <string>

// root header files
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TProfile.h"

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief plots ToT against spatial position of various hit pixels in a track-matched cluster
 *
 *
 */
class QShare1D: public TbAnalysis {

private:
	TH1D* h_tracks;
	TH1D* h_allHitsCenter;
	TH1D* h_nHits1Center;
	TH1D* h_nHits2Center;
	TH1D* h_nHits3PlusCenter;
	TH1D* h_allHitsNotCenter;
	TH2D* h_sumtot_allHitsCenter;
	TH2D* h_sumtot_nHits1Center;
	TH2D* h_sumtot_nHits2Center;
	TH2D* h_sumtot_nHits3PlusCenter;
	TH2D* h_tot_allHitsCenter;
	TH2D* h_tot_nHits1Center;
	TH2D* h_tot_nHits2Center;
	TH2D* h_tot_nHits3PlusCenter;
	TH2D* h_fractot_nHits2Center;

	void drawTotSharing(const TbConfig& config, const char* analysisName,
			const char* histoName, const char* drawOpts, TH2D* hsum,
			TH2D* hother, TH2D* hfrac) const;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //QSHARE1D_H
