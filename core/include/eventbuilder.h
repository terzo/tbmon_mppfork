#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H

class EventBuilder;
#include "event.h"
#include "tbconfig.h"

#include <string>

/*! \brief Mostly virtual base class for event builder processors
 *
 */
class EventBuilder {
public:
	virtual void init(TbConfig& config) {
	}
	;
	virtual void initRun(TbConfig &config) {
		;
	}
	;
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

	EventBuilder() :
			name("notset") {

	}
	;
	string name;
};

#endif
