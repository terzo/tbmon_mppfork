#include "batetacutter.h"

const int BatEtaCutter::batIdens[3] = { 81, 83, 86 };

void BatEtaCutter::init(TbConfig& config) {
	if (acceptRangeStart < 0.0 || acceptRangeStop > 1.0
			|| acceptRangeStart > acceptRangeStop) {
		cout << "[ etaCutter ]; Error: Invalid cut ranges" << endl;
		exit(-1);
	}

}

void BatEtaCutter::buildEvent(Event &event, map<int, Event>&, TbConfig &) {
	if (event.track == NULL) {
		cout << "[ etaCutter ]; Error: Enabled etaCutter for non-BAT data?!?"
				<< endl;
		exit(-1);
	}

	//Per-layer badness limits
	double badness[2];
	badness[0] = badness[1] = 0;

	//Loop over BAT planes
	for (int i = 0; i < numBats; i++) {
		int iden = batIdens[i];
		//Find cluster(s)
		for (int iCluster = 0; iCluster < event.track->nBatCluster;
				iCluster++) {
			BatCluster* cluster =
					(BatCluster*) event.track->batCluster[iCluster];
			if (cluster->iden == iden) {
				double eta = fmod((double) cluster->cstrp, 1.0);
				if (eta < acceptRangeStart) {
					badness[cluster->claye] += (acceptRangeStart - eta)
							* (acceptRangeStart - eta);
				} else if (eta > acceptRangeStop) {
					badness[cluster->claye] += (eta - acceptRangeStop)
							* (eta - acceptRangeStop);
				}
			}
		}
	}

	if (badness[0] > badness_limit[0] || badness[1] > badness_limit[1]) {
		event.fEtaCut = event::kBad;
		event.fTrack = event::kBad;
	} else
		event.fEtaCut == event::kGood;
}
