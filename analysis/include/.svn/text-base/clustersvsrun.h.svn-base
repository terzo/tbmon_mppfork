#ifndef CLUSTERSVSRUN_H
#define CLUSTERSVSRUN_H

// standard header files
#include <map>

// root header files
#include <TGraphErrors.h>

// tbmon header files
#include "event.h"
#include "tbanalysis.h"

/*! \brief calculates and plots mean cluster size (total, in x and in y) for each run
 *
 */
class ClustersVsRun: public TbAnalysis {

private:
	/// histograms shows mean total cluster size for each run
	TH1D* histo_matchClusterSizeVsRun;
	/// histograms shows mean total cluster in x size for each run
	TH1D* histo_matchClusterSizeXVsRun;
	/// histograms shows mean total cluster size in y for each run
	TH1D* histo_matchClusterSizeYVsRun;

	///total number of clusters per run
	std::map<int, int> nclusters;
	///sum of all cluster sizes per run
	std::map<int, int> sumclusters;
	///sum of all cluster sizes in x per run
	std::map<int, int> sumclusters_x;
	///sum of all cluster sizes in y per run
	std::map<int, int> sumclusters_y;

	///maximum cluster size on output plots
	double ClusterMax;
	///minimum cluster size on output plots
	double ClusterMin;

	/// graph shows mean total cluster size for each run
	TGraphErrors* h_matchClusterSizeVsRun;
	/// graph shows mean total cluster in x size for each run
	TGraphErrors* h_matchClusterSizeXVsRun;
	/// graph shows mean total cluster size in y for each run	///
	TGraphErrors* h_matchClusterSizeYVsRun;

	void drawAndSave(const TbConfig& config, const char* analysisName,
			const char* histoName, TGraphErrors* gr);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //CLUSTERSVSRUN_H
