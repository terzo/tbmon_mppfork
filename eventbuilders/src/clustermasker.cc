#include "clustermasker.h"

void ClusterMasker::buildEvent(Event &event, map<int, Event>&,
		TbConfig &config) {

	const DUT* dut = config.getDut(event.dut->getDUTid());
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();

	vector < vector<PllHit*> > tmpevent;
	tmpevent.clear();

	event.fTrackRegion = event::kGood;

	for (vector<vector<PllHit*> >::iterator ii = event.clusters.begin();
			ii != event.clusters.end(); ii++) {

		clusterIsGood = true;
		matchedTrack = false;

		for (vector<PllHit*>::iterator jj = ii->begin(); jj != ii->end();
				jj++) {
			int col((*jj)->col), row((*jj)->row);

			// check the cluster
			if (clusterIsGood) {

				pixelType = dut->maskedpixels[col + 1][row + 1];

				/*
				 cout << dut->getDUTid() << " " << col << " " << row << " " << pixelType  << endl;
				 cout << dut->maskedpixels[col+0][row+0] << " " << dut->maskedpixels[col+1][row+0] << " " << dut->maskedpixels[col+2][row+0] << endl;
				 cout << dut->maskedpixels[col+0][row+1] << " " << dut->maskedpixels[col+1][row+1] << " " << dut->maskedpixels[col+2][row+1] << endl;
				 cout << dut->maskedpixels[col+0][row+2] << " " << dut->maskedpixels[col+1][row+2] << " " << dut->maskedpixels[col+2][row+2] << endl;
				 cout <<  "--------" << endl;
				 */
				if (pixelType > 0) {
					clusterIsGood = false;
				}
			}

			// check the track
			if (fabs(event.trackX - (pitchX * col)) < 1.5 * pitchX
					&& fabs(event.trackY - (pitchY * row)) < 1.5 * pitchY) {
				matchedTrack = true;
			}
		}

		if (clusterIsGood) {
			tmpevent.push_back((*ii));
		} else if (matchedTrack && event.fTrackRegion == event::kGood) {
			event.fTrackRegion = event::kBad;
		}
	}

	event.clusters = tmpevent;
}

