#ifndef MAXCELLRESIDUALS_H
#define MAXCELLRESIDUALS_H

// root header files
#include <TH2D.h>

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Basic data quality checks; plot cut and un-cut X and Y residual, lvl1, track chi2, sensor and
 max ToT, etc. distributions; get approximate overall sensor efficiency numbers; this is a grab-bag, courtesy
 of Håvard
 *
 */
class MaxCellResiduals: public TbAnalysis {
private:
	TH1D* h_chi2;
	TH1D* h_lv1;
	TH1D* h_resX;
	TH1D* h_resY;
	TH1D* h_hitX;
	TH1D* h_hitY;
	TH1D* h_resXU;
	TH1D* h_resYU;
	TH1D* h_sumTot;
	TH1D* h_sensorTot;
	TH1D* h_matches;
	TH2D* h_hits;
	int m_tracks;
	int m_matches;
public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //MAXCELLRESIDUALS_H
