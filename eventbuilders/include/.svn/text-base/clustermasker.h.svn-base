#ifndef CLUSTERMAKER_H
#define CLUSTERMAKER_H

// standard header files
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"
#include "tbutils.h"

/*! \brief Goes through each cluster of an event (prerequisite: clusters have been built) and
 * checks for each pixel in a cluster if it is masked (either it is noisy or dead). In that case
 * the whole cluster is removed from the event and furthermore if the current track is matched the
 * trackregion (basically the current track) is marked as bad.
 *
 */
class ClusterMasker: public EventBuilder {
private:
	int pixelType; // Type of the pixel: 0 - good pixel, 1-9 dead pixel, 10 noisy pixel, <10 dead and noisy
//  bool clusterIsDead, clusterIsNoisy;
	bool clusterIsGood;
	bool matchedTrack;

public:

	ClusterMasker() {
		name = "ClusterMasker";

	}
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
};

#endif //CLUSTERMAKER_Hs
