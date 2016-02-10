#ifndef FULL3D_VADIM_H
#define FULL3D_VADIM_H

// standard header files
#include <math.h>

// root header files
#include "TGraph.h"

// tbmon header files
#include "simdut.h"
#include "event.h"
#include "simDutEdep.h"

/*! \brief This is a model for Full3D sensors, earlier presented by Vadim Kostioukhine
 *
 * This model treats the electrodes as cylinders
 * traversing the bulk perpendicularly,with a
 * efficiency fall-off
 *
 *    e(r) = 1 / (1 + exp( (R - r) / s))
 *
 * where R is the electrode radius, and s the
 * characteristic length of the fall-off.
 * *
 * TOT response proportional to edep.
 *
 * No timing response: always LVL1 = 1.
 */
class Full3D_Vadim: public Simdut {
private:
	double** qMatrix;    //Non-sparse storage of charge collected in pixels

	//List containing readout and bias electrodes x-position,
	// in "xmod" coordinates (see digitize())
	int N_readout, N_bias;
	double *readoutE_xpos, *biasE_xpos;

	/*
	 * Model parameters
	 */
	// TOT response
	/*
	 double totCalib;        //Conversion factor from MeV -> TOT
	 double chargeCalib;     //MeV pr. eh-pair
	 int    chargeThreshold; //Threhsold of preamp in #eh-pairs
	 double totThreshold;    // Same as above, but in TOT units
	 //   (calculated from above)
	 */
	double chargeThreshold; //Threshold of preamp in #eh-pairs
	double tot_a0, tot_a1, tot_a2;
	double tot_qMin, tot_qMax;
	double si_w;

	// Electrode parameters
	int nElec;        // Number of readout electrodes pr. pixel
	double R_readout; // Radius of electrodes
	double R_bias;
	double s_readout; // Characteristic fall-off of electrode efficiency
	double s_bias;

	//Squared distance between elecPos and (xmod,ymod)
	double elecDist2(double elecPosX, double elecPosY, double xmod,
			double ymod);

	//TOT response function
	inline int tot_response(double Q) {
		return (int) (tot_a0 * (tot_a1 + Q) / (tot_a2 + Q));
	}

public:
	Full3D_Vadim(const Event &event);
	virtual void init(TbConfig &config, const Event &event);
	virtual void digitize(Event& event, const TbConfig& config,
			std::vector<PllHit*>& hits);

	virtual void finalize(const TbConfig &config);

	~Full3D_Vadim();

};

#endif //FULL3D_VADIM_H
