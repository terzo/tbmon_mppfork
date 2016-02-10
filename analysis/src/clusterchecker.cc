#include "clusterchecker.h"

void ClusterChecker::init(TbConfig &config) {

	h_clusterMult = new TH1D("", ";Cluster Multiplicity", 20, 0, 20);
	h_nHits = new TH1D("", ";Number of Hits", 20, 0, 20);
	h_matchClusterSize = new TH1D("", ";Matched Cluster Size", 20, 0, 20);
	h_matchClusterSizeX = new TH1D("", ";Matched Cluster Size X", 20, 0, 20);
	h_matchClusterSizeY = new TH1D("", ";Matched Cluster Size Y", 20, 0, 20);
	h_unmatchClusterSize = new TH1D("", ";Un-Matched Cluster Size", 20, 0, 20);

}

void ClusterChecker::event(const TbConfig &config, const Event &event) {

	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		return;
	}
	// Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		return;
	}
	// Are we in a good region of the chip?
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);
	//if( matched == -1 ) { return; }

	int nclusters = event.clusters.size();
	h_clusterMult->Fill(nclusters);
	h_nHits->Fill(event.hits.size());

	int MaxSizeX;
	int MinSizeX;
	int MaxSizeY;
	int MinSizeY;

	for (int ii = 0; ii < event.clusters.size(); ii++) {
		if (ii == matched) {
			MinSizeX = 10000;
			MaxSizeX = 0;
			MinSizeY = 10000;
			MaxSizeY = 0;

			for (vector<PllHit*>::const_iterator jj =
					event.clusters.at(matched).begin();
					jj != event.clusters.at(matched).end(); ++jj) {
				if ((*jj)->col > MaxSizeX) {
					MaxSizeX = (*jj)->col;
				}
				if ((*jj)->col < MinSizeX) {
					MinSizeX = (*jj)->col;
				}
				if ((*jj)->row > MaxSizeY) {
					MaxSizeY = (*jj)->row;
				}
				if ((*jj)->row < MinSizeY) {
					MinSizeY = (*jj)->row;
				}
			}

			h_matchClusterSize->Fill(event.clusters.at(matched).size());
			h_matchClusterSizeX->Fill(MaxSizeX - MinSizeX + 1);
			h_matchClusterSizeY->Fill(MaxSizeY - MinSizeY + 1);

		} else {
			h_unmatchClusterSize->Fill(event.clusters.at(ii).size());
		}
	}

}

void ClusterChecker::finalize(const TbConfig &config) {

	config.drawAndSave(name, "multiplicity", h_clusterMult);
	config.drawAndSave(name, "pixelHits", h_nHits);
	config.drawAndSave(name, "matchClusterSize", h_matchClusterSize);
	config.drawAndSave(name, "matchClusterSizeX", h_matchClusterSizeX);
	config.drawAndSave(name, "matchClusterSizeY", h_matchClusterSizeY);
	config.drawAndSave(name, "unmatchClusterSize", h_unmatchClusterSize);

}
