//standard
#include <assert.h>

//root
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"

//tbmon
#include "clusters.h"
#include "clustersvsrun.h"

void ClustersVsRun::init(TbConfig & config) {

	TBALOG(kDEBUG3) << "Initializing ClustersVsRun processor." << std::endl;

	// get options or set standard values
	ClusterMin =
			config.cmdLineExtras_argGetter(
					(string) "A_clustersvsrun_MinSumTotClus", (double) 0.5,
					(string) "Minimum y-axis (mean cluster size) - may be overridden by automatic.");
	TBALOG(kDEBUG3) << "ClusterMin=" << ClusterMin << " has been set."
			<< std::endl;
	ClusterMax =
			config.cmdLineExtras_argGetter(
					(string) "A_clustersvsrun_MaxSumTotClus", (double) 1.5,
					(string) "Maximum y-axis (mean cluster size) - may be overridden by automatic.");
	TBALOG(kDEBUG3) << "ClusterMax=" << ClusterMax << " has been set."
			<< std::endl;

	for (std::vector<int>::const_iterator currentRun = config.runList.begin();
			currentRun != config.runList.end(); currentRun++) {

		nclusters[*currentRun] = 0;
		sumclusters[*currentRun] = 0;
		sumclusters_x[*currentRun] = 0;
		sumclusters_y[*currentRun] = 0;
		TBALOG(kDEBUG3) << "Initialized storage space for run "
				<< (int) (*currentRun) << "." << std::endl;
	};

	if (config.runList.size() < 2) {
		TBALOG(kERROR) << "Does not work properly with less than two runs."
				<< std::endl;
	};

	TBALOG(kDEBUG3)
			<< "Initialization of ClustersVsRun processor has been finished."
			<< std::endl;

}

void ClustersVsRun::event(const TbConfig &config, const Event &event) {

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

	int MaxSizeX;
	int MinSizeX;
	int MaxSizeY;
	int MinSizeY;

	//loop over all clusters in this event
	for (int ii = 0; ii < event.clusters.size(); ii++) {
		TBALOG(kDEBUG3) << "ClustersVsRun looks at cluster " << ii
				<< " of event " << config.currentEntry << "." << std::endl;
		//only look at matched cluster
		if (ii == matched) {
			TBALOG(kDEBUG3) << "ClustersVsRun looks at matched cluster " << ii
					<< " of event " << config.currentEntry << "." << std::endl;
			MinSizeX = 10000;
			MaxSizeX = 0;
			MinSizeY = 10000;
			MaxSizeY = 0;

			//loop over all hits in this cluster
			for (vector<PllHit*>::const_iterator jj =
					event.clusters.at(matched).begin();
					jj != event.clusters.at(matched).end(); ++jj) {
				TBALOG(kDEBUG3) << "ClustersVsRun looks at hit "
						<< std::distance(event.clusters.at(matched).begin(), jj)
						<< " in cluster " << ii << " of event "
						<< config.currentEntry << "." << std::endl;
				if ((*jj)->col > MaxSizeX) {
					MaxSizeX = (*jj)->col;
				}
				TBALOG(kDEBUG3) << "ClustersVsRun: MaxSizeX= " << MaxSizeX
						<< "." << std::endl;
				if ((*jj)->col < MinSizeX) {
					MinSizeX = (*jj)->col;
				}
				TBALOG(kDEBUG3) << "ClustersVsRun: MinSizeX= " << MinSizeX
						<< "." << std::endl;
				if ((*jj)->row > MaxSizeY) {
					MaxSizeY = (*jj)->row;
				}
				TBALOG(kDEBUG3) << "ClustersVsRun: MaxSizeY= " << MaxSizeY
						<< "." << std::endl;
				if ((*jj)->row < MinSizeY) {
					MinSizeY = (*jj)->row;
				}
				TBALOG(kDEBUG3) << "ClustersVsRun: MinSizeY= " << MinSizeY
						<< "." << std::endl;
			}

			int cluster = event.clusters.at(matched).size();
			TBALOG(kDEBUG3) << "ClustersVsRun: Global cluster size= "
					<< cluster << "." << std::endl;
			nclusters[config.currentRun]++;
			TBALOG(kDEBUG3) << "ClustersVsRun: Number of clusters= "
					<< nclusters[config.currentRun] << "." << std::endl;
			sumclusters[config.currentRun] = sumclusters[config.currentRun]
					+ cluster;
			TBALOG(kDEBUG3) << "ClustersVsRun: All clusters added up= "
					<< sumclusters[config.currentRun] << "." << std::endl;
			sumclusters_x[config.currentRun] = sumclusters_x[config.currentRun]
					+ (MaxSizeX - MinSizeX + 1);
			TBALOG(kDEBUG3) << "ClustersVsRun: MaxSizeX - MinSizeX + 1= "
					<< (MaxSizeX - MinSizeX + 1) << "." << std::endl;
			sumclusters_y[config.currentRun] = sumclusters_y[config.currentRun]
					+ (MaxSizeY - MinSizeY + 1);
			TBALOG(kDEBUG3) << "ClustersVsRun: MaxSizeY - MinSizeY + 1= "
					<< (MaxSizeY - MinSizeY + 1) << "." << std::endl;

		} else {
			//Can add option to output unmatched cluster sizes.
		}
	}

}

void ClustersVsRun::finalize(const TbConfig &config) {

	int nRuns = nclusters.size();

	double* runs = new double[nRuns];
	double* clustermatched = new double[nRuns];
	double* clustermatched_x = new double[nRuns];
	double* clustermatched_y = new double[nRuns];

	histo_matchClusterSizeVsRun = new TH1D("",
			";Run Number;(Matched) Cluster Size",
			config.runList.at(nRuns - 1) - config.runList.at(0),
			config.runList.at(0), config.runList.at(nRuns - 1));
	histo_matchClusterSizeXVsRun = new TH1D("",
			";Run Number;(Matched) Cluster Size X",
			config.runList.at(nRuns - 1) - config.runList.at(0),
			config.runList.at(0), config.runList.at(nRuns - 1));
	histo_matchClusterSizeYVsRun = new TH1D("",
			";Run Number;(Matched) Cluster Size Y",
			config.runList.at(nRuns - 1) - config.runList.at(0),
			config.runList.at(0), config.runList.at(nRuns - 1));

	for (int irun = 0; irun < nRuns; irun++) {
		int run = config.runList.at(irun);
		runs[irun] = (double) run;

		if (nclusters[run] != 0 && sumclusters[run] != 0) {
			// output of cluster sizes as text
			TBALOG(kINFO) << "Clusters for run " << run << ": " << endl;
			TBALOG(kINFO) << "Mean (matched) cluster size: "
					<< (double) sumclusters[run] / nclusters[run] << endl;
			TBALOG(kINFO) << "Mean (matched) cluster x size: "
					<< (double) sumclusters_x[run] / nclusters[run] << endl;
			TBALOG(kINFO) << "Mean (matched) cluster y size: "
					<< (double) sumclusters_y[run] / nclusters[run] << endl;

			clustermatched[irun] = (double) sumclusters[run] / nclusters[run];// Save mean of matched cluster size
			clustermatched_x[irun] = (double) sumclusters_x[run]
					/ nclusters[run];// Save mean of matched cluster size in x direction
			clustermatched_y[irun] = (double) sumclusters_y[run]
					/ nclusters[run];// Save mean of matched cluster size in y direction

			if (ClusterMin > clustermatched[irun])
				ClusterMax = clustermatched[irun] - 1; // Change plot axis max
			if (ClusterMax < clustermatched[irun])
				ClusterMax = clustermatched[irun] + 1;	// Change plot axis min

		} else {
			cout << "No clusters found." << endl;
			clustermatched[irun] = 0;
		}

		histo_matchClusterSizeVsRun->SetBinContent(
				run - config.runList.at(0) + 1, clustermatched[irun]);
		histo_matchClusterSizeXVsRun->SetBinContent(
				run - config.runList.at(0) + 1, clustermatched_x[irun]);
		histo_matchClusterSizeYVsRun->SetBinContent(
				run - config.runList.at(0) + 1, clustermatched_y[irun]);

	}

	// Save histograms to root file and to folder

	config.drawAndSave(this->name, (const char*) "matchClusterSizeVsRun",
			histo_matchClusterSizeVsRun);
	config.drawAndSave(this->name, (const char*) "matchClusterSizeXVsRun",
			histo_matchClusterSizeXVsRun);
	config.drawAndSave(this->name, (const char*) "matchClusterSizeYVsRun",
			histo_matchClusterSizeYVsRun);

	h_matchClusterSizeVsRun = new TGraphErrors(nRuns, runs, clustermatched, 0,
			0);
	h_matchClusterSizeXVsRun = new TGraphErrors(nRuns, runs, clustermatched_x,
			0, 0);
	h_matchClusterSizeYVsRun = new TGraphErrors(nRuns, runs, clustermatched_y,
			0, 0);

	drawAndSave(config, name, "matchClusterSizeVsRun", h_matchClusterSizeVsRun);
	drawAndSave(config, name, "matchClusterSizeXVsRun",
			h_matchClusterSizeXVsRun);
	drawAndSave(config, name, "matchClusterSizeYVsRun",
			h_matchClusterSizeYVsRun);

	delete[] clustermatched;

}

void ClustersVsRun::drawAndSave(const TbConfig& config,
		const char* analysisName, const char* histoName, TGraphErrors* gr) {
	std::string tempExtension;
	tempExtension = config.plotExtension;
	if (tempExtension != ".none") {
		char* fileName = config.buildHistName(analysisName, histoName);
		const DUT* dut = config.getDut(this->iden);

		TCanvas* c_clust = new TCanvas("c_clust", "c_clust");
		c_clust->cd();

		gr->Draw("AP");
		gr->GetXaxis()->SetTitle("Run number");
		gr->GetYaxis()->SetTitle("Matched Cluster Size");
		gr->GetYaxis()->SetRangeUser(ClusterMin, ClusterMax);
		c_clust->SaveAs(fileName);

		delete c_clust;
	};
}
