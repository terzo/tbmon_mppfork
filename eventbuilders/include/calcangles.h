#ifndef CALCANGLES_H
#define CALCANGLES_H

// root header files
#include <string>

// tbmon header files
#include "eventbuilder.h"

/*! \brief Class calculating the actual rotation of the DUT using information from
 *  the GEAR file and the track fitting during reconstruction
 */
class CalcAngles: public EventBuilder {

public:
	virtual void init(TbConfig& config) {
	}
	;
	virtual void initRun(TbConfig &config);
	virtual void initEvent(TbConfig &config, map<int, Event> &events) {
		;
	}
	;
	virtual void buildEvent(Event &event, map<int, Event> &events,
			TbConfig &config) {
	}
	;
	virtual void finalizeRun(TbConfig &config) {
		;
	}
	;
	virtual void finalize(TbConfig &config) {
		;
	}
	;

	CalcAngles() {
		name = "CalcAngles";
	}
	;
};

#endif //CALCANGLES_H
