#ifndef LVL1CUTS_H
#define LVL1CUTS_H

// standard header files
#include <cmath>
#include <vector>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief transfers hits from rawHits list to hits if they fulfill the lv1 cut criterium defined in DUT object
 *
 */
class LVL1Cuts: public EventBuilder {
public:
	LVL1Cuts() {
		name = "LVL1Cuts";
	}
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
};

#endif //LVL1CUTS_H
