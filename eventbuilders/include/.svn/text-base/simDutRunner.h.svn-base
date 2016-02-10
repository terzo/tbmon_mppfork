#ifndef SIMDUTRUNNER_H
#define SIMDUTRUNNER_H

#include <iostream>
#include <map>
#include <stdlib.h>

#include "eventbuilder.h"
#include "simdut.h"

//DUT models
#include "pixel_simple.h"
#include "Full3D_HP.h"
#include "Full3D_Vadim.h"

/*!  \brief SIMULATION: ?
 *
 */
class simDutRunner: public EventBuilder {
private:
	string modelName;
	Simdut* dutSim;
public:

	virtual void init(TbConfig& config, const Event &event);
	virtual void buildEvent(Event &event, map<int, Event> &events,
			TbConfig &config);
	virtual void finalize(TbConfig& config);
};

#endif //SIMDUTRUNNER_H
