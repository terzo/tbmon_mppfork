#ifndef SIMRESIDUALS_H
#define SIMRESIDUALS_H


#include "tbanalysis.h"
#include "tbconfig.h"

#include <TH1D.h>
#include <TH2D.h>

/*! \brief SIMULATION: Very simple analysis comparing truth information with track extrapolated into the device plane.
 *
 * Yields track resolution in DUT.
 *
 * Kyrre Sjobak
 */
class simResiduals: public TbAnalysis {
private:
	TH1D* h_resX;
	TH1D* h_resY;

	//Truth vs track correlations
	TH2D* h_XvsX;
	TH2D* h_XvsY;
	TH2D* h_YvsX;
	TH2D* h_YvsY;

public:

	//Stuff to be done before anything else
	virtual void init(TbConfig &config);
	//How to process an event
	virtual void event(const TbConfig &config, const Event &event);
	//Stuff to be done after all events have been processed
	virtual void finalize(const TbConfig &config);
	//Stuff to be done at the beginning of each run
	//virtual void initRun(const TbConfig &config){;}
	//Stuff to be done at the end of each run
	//virtual void finalizeRun(const TbConfig &config){;}
};

#endif //SIMRESIDUALS_H
