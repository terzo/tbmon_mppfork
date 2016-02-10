#ifndef CLUSTERCHECKER_H
#define CLUSTERCHECKER_H

// root header files
#include <TH1D.h>

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"
#include "tbconfig.h"

/*! \brief basic cluster analysis - plots track-matched and un-matched cluster size (number of cells),
cluster multiplicity and number of pixel hits per event.
 *
 */
class ClusterChecker: public TbAnalysis {
private:
	TH1D* h_clusterMult;
	TH1D* h_nHits;
	TH1D* h_matchClusterSize;
	TH1D* h_matchClusterSizeX;
	TH1D* h_matchClusterSizeY;
	TH1D* h_unmatchClusterSize;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //CLUSTERCHECKER_H
