#ifndef CHECKCENTRALREGION_H
#define CHECKCENTRALREGION_H

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

/*! \brief Checks whether the track is in the central region of the sensor and marks the event otherwise as bad
 * (TODO: hardcoded values, this shouldn't probably be done this way as any processor wanting to analyze stuff in the edge
 * won't work, depends on addMaskedPixels() --> two competing ways of masking!)
 *
 */
class CheckCentralRegion: public EventBuilder {
public:
	CheckCentralRegion() {
		name = "CheckCentralRegion";
	}
	virtual void initRun(TbConfig &config) {
		;
	}
	;
	virtual void buildEvent(Event &event, map<int, Event> &events, TbConfig &);
};

#endif //CHECKCENTRALREGION_H
