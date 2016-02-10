#ifndef MASKANDLVL1_H
#define MASKANDLVL1_H

// standard header files
#include <cmath>
#include <vector>
#include <stdlib.h>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief transfers hits from rawHits list to hits, according to the mask of the DUT, and performs a command line lvl1 cut
 *
 *  syntax for command line lvl1 input: -P:B_maskandlvl1_lvl1values "5 6 8" (for example)
 */
class MaskAndLvl1: public EventBuilder {

public:
	MaskAndLvl1() {
		name = "MaskAndLvl1";
	}
	virtual void init(TbConfig &);
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);

private:
	vector<int> lvl1;
	bool lvl1cut;
};

#endif //MASKANDLVL1_H
