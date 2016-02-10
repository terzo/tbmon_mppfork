#ifndef CHECKALIGN_H
#define CHECKALIGN_H

// standard header files
#include <vector>

// root header files
#include <TProfile.h>
#include <TH1D.h>

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Plots average residual vs hit position (x vs y, and all combinations) to assess whether the
sensors were well-aligned during track reconstruction; also quickly spot runs that don’t appear to be well
aligned by checking their mean x and y residual values; prints out suspicious runs. Generates a (x,y) vector per run
which is written into a file and can be read back by the translation event builder to shift each track point.
 *
 */
class CheckAlign: public TbAnalysis {

private:
	bool doCuts;
	TProfile* h_hitxVresx;
	TProfile* h_hityVresx;
	TProfile* h_hitxVresy;
	TProfile* h_hityVresy;
	TH1D* h_resX;
	TH1D* h_resY;
	//DEBUG HISTOGRAM
	TH1D* h_cuts;
	vector<int> badRuns;
	vector<double> deltaX;
	vector<double> deltaY;
	vector<int> run;

public:
	virtual void init(TbConfig &config);
	virtual void initRun(const TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalizeRun(const TbConfig &config);
	virtual void finalize(const TbConfig &config);
};

#endif //CHECKALIGN_H
