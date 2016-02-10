#ifndef CLUSTERFINDER_H
#define CLUSTERFINDER_H

// standard header files
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief Event builder which takes rawHits and Hits and creates rawClusters and Clusters.
 *
 */
class ClusterFinder: public EventBuilder {

private:
	//
	int addNeighbors(vector<PllHit*> &cluster, list<PllHit*> &hits);

public:
	//
	ClusterFinder() {
		name = "ClusterFinder";
	}
	//
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
	//
	void findClusters(vector<PllHit*> &hits,
			vector<vector<PllHit*> > &clusters);
};

#endif //CLUSTERFINDER_H
