#ifndef CLUSTERPIXTOTDISTR_H
#define CLUSTERPIXTOTDISTR_H

// standard header files
#include <string>
#include <map>

// tbmon header files
#include "tbanalysis.h"

class TTree;
class TFile;
class TH1F;

/*! \brief Saves and plots ToT vs. position of pixel in cluster and 1st pixel ToT
 *
 */
class ClusterPixToT: public TbAnalysis {
private:
	/// Tree to save ToT information of pixels
	TTree *_outTree;

	/// Stores ToT information [minPosToT][_pixelPosX][_pixelToT]
	std::map<int, std::map<int, std::map<int, int> > > _1stPixToTdependentToT;
	/// Stores ToT distribution of single pixels in cluster
	std::map<int, TH1F*> _singlePixToTDistr;
	/// Stores ToT distribution dependent on 1st pixel ToT and pixel position
	std::map<int, TH1F*> _1stPixDepTotDistr;
	/// Stores ToT distribution of single pixel dependent on 1st pixel ToT
	std::map<int, TH1F*> _1stPixDepSinglePixTotDistr;
	///
	int _clusterNumber;
	/// Calculated cluster size in X
	int _clusterSizeX;
	/// Calculated cluster size in Y
	int _clusterSizeY;
	/// Pixel ToT
	int _pixelToT;
	/// X position of pixel in cluster
	int _pixelPosX;
	/// Y position of pixel in cluster
	int _pixelPosY;

	/// Only consider cluster with a certain X size for plots
	int _fixedClusterSizeX;
	/// Only consider cluster with a certain Y size for plots. Currently fixed to 1
	int _fixedClusterSizeY;

public:
	ClusterPixToT() :
			_outTree(0), _fixedClusterSizeY(1) {
	}
	;
	virtual ~ClusterPixToT() {
	}
	;
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config) {
	}
	;
	virtual void initRun(const TbConfig &config);
	virtual void finalizeRun(const TbConfig &config);
};

#endif //CLUSTERPIXTOTDISTR_H
