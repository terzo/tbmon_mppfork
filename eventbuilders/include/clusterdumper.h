#ifndef CLUSTERDUMPER_H
#define CLUSTERDUMPER_H

// standard header files
#include <vector>
#include "stdlib.h"

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*!  \brief Event builder which reads a list of cluster sizes from the command line. Then it goes
 * through each already built cluster and removes all that are not of one of the sizes previously
 * defined.
 *
 */
class ClusterDumper: public EventBuilder {

public:
	ClusterDumper() {
		name = "ClusterDumper";
		dump = false;
	}
	virtual void init(TbConfig &config);
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
	void findClusters(vector<PllHit*> &hits,
			vector<vector<PllHit*> > &clusters);

private:
	bool dump;
	vector<int> sizes;
	bool checkSize(int clusize);
};

#endif //CLUSTERDUMPER_H
