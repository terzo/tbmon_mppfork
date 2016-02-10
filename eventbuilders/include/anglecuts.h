#ifndef ANGLECUTS_H
#define ANGLECUTS_H

// standard header files
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief Applies angle cuts where the paramters are read from a file
 *
 */
class AngleCuts: public EventBuilder {
private:
	/*! \brief angle data for one run
	 *
	 */
	struct AngleStats {
		int run;
		double meanX;
		double sigmaX;
		double meanY;
		double sigmaY;
	};
	map<int, AngleStats> angleStats;
	AngleStats currentCuts;
	double sigmas;
public:
	AngleCuts(const char* fileName = NULL, double sigmas = 1.5);
	virtual void initRun(TbConfig &config);
	virtual void buildEvent(Event &event, map<int, Event> &events, TbConfig &);
};

#endif //ANGLECUTS_H
