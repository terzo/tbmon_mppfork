#ifndef CORRELATIONS_H
#define CORRELATIONS_H

// standard header files
#include <string>

// root header files
#include "TH2D.h"

// tbmon header files
#include "event.h"
#include "tbanalysis.h"

/*! \brief Plots correlations of hits with tracks
 *
 */
class Correlation: public TbAnalysis {
private:
	TH2D* h_colvx;
	TH2D* h_colvy;
	TH2D* h_rowvx;
	TH2D* h_rowvy;

	TH2D* h_occup;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //CORRELATIONS_H
