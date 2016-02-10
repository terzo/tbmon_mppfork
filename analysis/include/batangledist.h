#ifndef BATANGLEDIST_H
#define BATANGLEDIST_H

// standard header files
#include <cmath>
#include <fstream>
#include <iostream>

// root header files
#include <TH1D.h>

// tbmon header files
#include "event.h"
#include "tbanalysis.h"
#include "Track.h"

/*! \brief ONLY FOR BAT DATA (no longer maintained - will be probably removed at some point): plots angle distribution and consequences of angle cut
 *
 */
class BatAngleDist: public TbAnalysis {
private:
	TH1D* h_angleXTmp;
	TH1D* h_angleYTmp;
	TH1D* h_acceptedAnglesX;
	TH1D* h_rejectedAnglesX;
	TH1D* h_acceptedAnglesY;
	TH1D* h_rejectedAnglesY;
	ofstream stream;
public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
	virtual void initRun(const TbConfig &config);
	virtual void finalizeRun(const TbConfig &config);
};

#endif //BATANGLEDIST_H
