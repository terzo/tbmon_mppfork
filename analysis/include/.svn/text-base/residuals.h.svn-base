#ifndef RESIDUALS_H
#define RESIDUALS_H

#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

//standard
#include <algorithm>
#include <cmath>
#include <map>
#include <string>

//root
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveStats.h"

/*! \brief This class creates so called residual plots which are useful to check the quality of the reconstruction and analysis as well as the telescope resolution. Basically a residual
 * is the difference between the calculated cluster center and the traversing point of its matched track.
 *
 */
class Residuals: public TbAnalysis {

private:

	double pitchX; ///< sensor pitch in central region in x-direction (copied from DUT)
	double pitchY; ///< sensor pitch in central region in y-direction (copied from DUT)

	TH1D* h_angleX; ///< Histogram of dx/dz from the eutracks tree. Only the first track in multi track event is handled and a cut is made on the nmatch criterium.
	TH1D* h_angleY; ///< Histogram of dy/dz from the eutracks tree. Only the first track in multi track event is handled and a cut is made on the nmatch criterium.
	TH1D* h_clusizeX; ///< Cluster size distribution of all clusters in x (columns) which passed the following cuts: central region and matched
	TH1D* h_clusizeY; ///< Cluster size distribution of all clusters in y (rows) which passed the following cuts: central region and matched
	TH1D* h_clusize; ///< Total cluster size distribution of all clusters which passed the following cuts: central region and matched
	TH2D* h_clusizeXY; ///< Cluster size distribution of all clusters in x-y which passed the following cuts: central region and matched
	TH1D* h_ecorrsX; ///< Distribution of eta correction in x
	TH1D* h_ecorrsY; ///< Distribution of eta correction in y

	map<int, vector<TH1D*> > hhh_res; ///< XY residual plots calculated with four different cluster position algorithm and a virtually unlimited number of cluster sizes
	map<int, vector<TH1D*> > hhh_resX; ///< X residual plots calculated with four different cluster position algorithm and a virtually unlimited number of cluster sizes
	map<int, vector<TH1D*> > hhh_resY; ///< Y residual plots calculated with four different cluster position algorithm and a virtually unlimited number of cluster sizes

	bool use_minimum_sum_of_tot_per_cluster; ///< Says whether a cut on the minimum ToT sum of a cluster shall be applied.
	int minimum_sum_of_tot_per_cluster; ///< If a cluster does not contain at least this number of entries, it is not used for the residual calculation.
	bool use_minimum_entries_to_plot_histo; ///< Says whether a cut on the minimum number of entries in a plot to be plotted shall be applied.
	int minimum_entries_to_plot_histo; ///< If a residual plot does not contain at least this number of entries, it is not drawn to a file and included in the root file during the finalize method.

	void fitResiduals(const TbConfig& config, const char* analysisName,
			const char* histoName, TH1 *h_residuals, double alfPitch,
			double xmin, double xmax);

	void fitAllClusResiduals(const TbConfig& config, const char* analysisName,
			const char* histoName, TH1 *h_residuals, double alfPitch,
			double xmin, double xmax);

	void makePlot(TH1 *h_residuals, char * nameString);
	void getResiduals(const vector<PllHit*> &cluster, const Event &event,
			vector<TH1D*> &histosX, vector<TH1D*> &histosY,
			vector<TH1D*> &histosAll);
	void compareResidualsForOneClustersize(const TbConfig& config,
			const char* analysisName, const char* histoName, vector<TH1D*>& hh);
	void checkEtaCorrections(const vector<PllHit*> &cluster,
			const Event &event);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //RESIDUALS_H
