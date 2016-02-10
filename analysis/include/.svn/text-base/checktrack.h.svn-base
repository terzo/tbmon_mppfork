#ifndef CHECKTRACK_H
#define CHECKTRACK_H

// standard header files
#include <vector>

// root header files
#include <TH2D.h>

// tbmon header files
#include "event.h"
#include "tbanalysis.h"
#include "tbutils.h"

/*! \brief Plots chi2 distribution and badRegion map.
 *
 */
class CheckTrack: public TbAnalysis {
private:
	TH1D* h_chi2;
	TH2D* h_badRegion;
public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //CHECKTRACK_H
