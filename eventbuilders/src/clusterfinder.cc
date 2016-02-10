#include "clusterfinder.h"

int ClusterFinder::addNeighbors(vector<PllHit*> &cluster, list<PllHit*> &hits) {
	int counter(0); //How many hits are added to the cluster this iteration?
	for (list<PllHit*>::iterator hit = hits.begin(); hit != hits.end(); hit++) {
		for (vector<PllHit*>::iterator it = cluster.begin();
				it != cluster.end(); it++) {
			//Everything within a 3x3 area around a hit is assumed to belong in the cluster
			if (fabs((*it)->row - (*hit)->row) > 1) {
				continue;
			}
			if (fabs((*it)->col - (*hit)->col) > 1) {
				continue;
			}
			cluster.push_back((*hit));
			hit = hits.erase(hit);
			counter++;
			break;
		}
	}
	return (counter);
}
void ClusterFinder::buildEvent(Event &event, map<int, Event> &, TbConfig &) {
	//Make a copy of all PllHits* to a list of available hits
	findClusters(event.hits, event.clusters);
	findClusters(event.rawHits, event.rawClusters);
	event.fClusters = event::kGood;
	event.fEtaCorrections = event::kGood; //Check if etacorrections are available
	if (event.dut->ecorrs.eCorrs.ecorrX.size() != 100) {
		event.fEtaCorrections = event::kBad;
	}
	if (event.dut->ecorrs.eCorrs.ecorrY.size() != 100) {
		event.fEtaCorrections = event::kBad;
	}
}

void ClusterFinder::findClusters(vector<PllHit*> &hits,
		vector<vector<PllHit*> > &clusters) {
	//Make a copy of all PllHits* to a list of available hits
	list<PllHit*> tmphits;
	tmphits.resize(hits.size());
	copy(hits.begin(), hits.end(), tmphits.begin());
	//Loop until all hits are put into a cluster
	while (not tmphits.empty()) {
		vector<PllHit*> cluster;
		//Cluster seed is first unused hit
		cluster.push_back(tmphits.front());
		tmphits.pop_front();
		// Add neighbours until cluster stops growing
		while (addNeighbors(cluster, tmphits) > 0) {
			;
		}
		clusters.push_back(cluster);
	}
}
