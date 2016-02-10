#ifndef BATTRACK_H
#define BATTRACK_H

// standard header files
#include <map>

// root header files
#include <TBranch.h>
#include "TTree.h"

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief ONLY FOR BAT DATA (no longer maintained - will be probably removed at some point): Reads in BAT data from file
 *
 */
class BatTrack: public EventBuilder {
private:
	TTree* tracktree;
	Track* m_track;
public:
	BatTrack() {
		name = "BatTrack";
	}
	;

	virtual void finalizeRun(TbConfig& config);
	virtual void initRun(TbConfig &config);
	virtual void buildEvent(Event &event, map<int, Event> &events, TbConfig &);
	virtual void initEvent(TbConfig &config, map<int, Event> &events);
	Track* getTrack() {
		return m_track;
	}
	;
};

#endif //BATTRACK_H
