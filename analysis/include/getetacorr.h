#ifndef GETETACORR_H
#define GETETACORR_H

// standard header files
#include <cmath>
#include <iostream>
#include <fstream>

// root header files
#include <TH1D.h>

// tbmon header files
#include <clusters.h>
#include "event.h"
#include "tbanalysis.h"

/*! \brief calculates eta correction for each DUT and run (valid for one and two hit clusters) and plots charge weighted residuals in x and y
 *
 */
class GetEtaCorr: public TbAnalysis {

private:
	TH1D* h_chargeWeightedX; ///< charge weighted cluster x position
	TH1D* h_chargeWeightedY; ///< charge weighted cluster y position
	TH1D* h_etacorrsX; ///< eta correction in x
	TH1D* h_etacorrsY; ///< eta correction in y

	virtual void printEtaCalib(TH1D* histo, TH1D* etacorrs, ofstream &stream);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //GETETACORR_H
