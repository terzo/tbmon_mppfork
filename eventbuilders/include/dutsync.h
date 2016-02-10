#ifndef DUTSYNCH_H
#define DUTSYNCH_H

// standard header files
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"

/*! \brief Event builder which checks for desynch? Might be for BAT data as some strange file is read in.
 *
 */
class DutSync: public EventBuilder {
private:
	map<int, vector<bool> > syncMap;
public:
	DutSync() {
		name = "DutSync";
	}
	bool checkCurrent;
	vector<bool> acceptSync;
	void readFile(const char* fileNeam);
	virtual void initRun(TbConfig &config);
	virtual void buildEvent(Event &event, map<int, Event> &events,
			TbConfig &config);
};

#endif //DUTSYNCH_H
