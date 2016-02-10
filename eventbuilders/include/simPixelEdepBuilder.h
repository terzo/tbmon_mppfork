#ifndef SIMPIXELEDEPBUILDER_H
#define SIMPIXELEDEPBUILDER_H

#include <fstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <string>

#include "event.h"
#include "eventbuilder.h"
#include "simPixelEdep.h"

/*! \brief SIMULATION: This eventbuilder reads edep data from TestBeamSim output and stores in in event.simData->
 *
 */
class simPixelEdepBuilder: public EventBuilder {
private:
	map<int, ifstream*> pixelFiles;
	//map<int,vector<simPixelEdep> > edeps;

	char* bitBucket;

public:
	simPixelEdepBuilder() {
		name = "simPixelEdepBuilder";
		bitBucket = new char[256];
	}
	~simPixelEdepBuilder() {
		delete bitBucket;
	}

	virtual void initRun(TbConfig &config);
	virtual void initEvent(TbConfig &config, map<int, Event> &events);
	virtual void buildEvent(Event &event, map<int, Event> &events,
			TbConfig& config);
	virtual void finalizeRun(TbConfig& config);
};

#endif //SIMPIXELEDEPBUILDER_H
