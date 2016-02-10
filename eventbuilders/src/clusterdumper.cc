#include "clusterdumper.h"

void ClusterDumper::init(TbConfig &config) {
	string cluststring = config.cmdLineExtras_argGetter("B_clusterdumper",
			string(""), "Cluster sizes");
	if (cluststring == "")
		dump = false;
	else {
		dump = true;
		string ints;
		stringstream cluststream(cluststring);
		while (getline(cluststream, ints, ' ')) {
			sizes.push_back(atoi(ints.c_str()));
		}
		//TBALOG (kDEBUG) << "Reading in lvl1 cut values: ";
		//for(vector<int>::const_iterator i = lvl1.begin(); i != lvl1.end(); i++) {TBALOG (kDEBUG) << (*i) << " ";}
	}

}

void ClusterDumper::buildEvent(Event &event, map<int, Event>&, TbConfig &) {
	if (!dump)
		return;
	vector<event::cluster_t> properClusters;
	for (vector<event::cluster_t>::const_iterator i = event.clusters.begin();
			i != event.clusters.end(); i++) {
		if (checkSize(i->size()))
			properClusters.push_back(*i);
	}
	event.clusters = properClusters;
}

bool ClusterDumper::checkSize(int clusize) {
	for (vector<int>::const_iterator i = sizes.begin(); i != sizes.end(); i++) {
		if (clusize == (*i))
			return true;
	}
	return false;
}
