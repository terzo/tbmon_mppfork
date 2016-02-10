#ifndef SIMDUTEDEP_H
#define SIMDUTEDEP_H

#include "tbanalysis.h"
#include "tbconfig.h"

#include <TH1D.h>

/*! \brief SIMULATION: Simple "test" analysis for checking that simulated edep info in DUTs work well.
 *
 * Kyrre Sjobak
 */
class simDutEdep: public TbAnalysis {
private:
	TH1D* h_sumEdep;

	//Wants: Measure spread of edeps =>
	//max achivable resolution
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

#endif //SIMDUTEDEP_H
