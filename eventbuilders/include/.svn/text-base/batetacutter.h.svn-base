#ifndef BATETACUTTER_H
#define BATETACUTTER_H

/*
 * This experimental builder allows for
 * setting track cuts, based on the value of
 * eta in the BAT planes.
 *
 * Experimental, may lead to strange
 * artifacts in beam profiles.
 *
 * Kyrre Ness Sjøbæk
 * k.n.sjobak@fys.uio.no
 */

// standard header files
#include <list>
#include <vector>
#include <sstream>
#include <stdlib.h>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief ONLY FOR BAT DATA (no longer maintained - will be probably removed at some point): Track cuts based on eta.
 *
 */
class BatEtaCutter: public EventBuilder {
private:
	double acceptRangeStart;
	double acceptRangeStop;

	static const int batIdens[3];
	static const int numBats = 3;

	double badness_limit[2];

public:
	//This builder gets its default value from constructor argument
	BatEtaCutter(double in_acceptRangeStart, double in_acceptRangeStop) :
			acceptRangeStart(in_acceptRangeStart), acceptRangeStop(
					in_acceptRangeStop) {
		name = "etaCutter";

		badness_limit[0] = badness_limit[1] = 0.05;
	}
	;
	virtual void init(TbConfig& config);
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
};

#endif //BATETACUTTER_H
