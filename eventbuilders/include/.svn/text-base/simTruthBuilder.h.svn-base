#ifndef SIMTRUTHBUILDER_H
#define SIMTRUTHBUILDER_H

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "event.h"
#include "eventbuilder.h"
#include "simTruthHit.h"

using namespace std;

/*! \brief SIMULATION: This eventbuilder reads truth data from TestBeamSim output and stores it in event.simData->truthHits[]
 *
 * Kyrre N. Sjobak
 */
class simTruthBuilder: public EventBuilder {
private:
	ifstream truthFile;

	vector<simTruthHit*> truthHits;

	char* bitBucket;

	bool useTruthTracking;

public:
	simTruthBuilder() {
		name = "simTruthBuilder";
		bitBucket = new char[256];
	}
	~simTruthBuilder() {
		delete bitBucket;
	}

	virtual void init(TbConfig &config);
	virtual void initRun(TbConfig &config);
	virtual void initEvent(TbConfig &config, const Event &event,
			map<int, Event> &events);
	virtual void buildEvent(Event &event, map<int, Event> &events, TbConfig &);
	virtual void finalizeRun(TbConfig& config);

};

#endif //SIMTRUTHBUILDER_H
