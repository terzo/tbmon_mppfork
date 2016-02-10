#ifndef FULL3D_HP_H
#define FULL3D_HP_H

// root header files
#include "TGraph.h"

// tbmon header files
#include "simdut.h"
#include "event.h"
#include "simDutEdep.h"

/*! \brief This is a simple model for Full3D sensors (HP = "HolePunch").
 *
 * It treats the electrodes as cylinders of radius R with
 * zero-efficiency, traversing the bulk perpendicularly.
 *
 * TOT response proportional to edep, always LVL1 = 1.
 */
class Full3D_HP: public Simdut {
private:
	double** qMatrix;    //Non-sparse storage of charge collected in pixels

	//List of 2-vectors containing readout and bias electrodes position,
	// in "xmod" coordinates (see digitize())
	double* readoutE_xpos;
	double* biasE_xpos;
	/*
	 * Model parameters
	 */
	// TOT response
	/* OLD linear model
	 double totCalib;        //Conversion factor from MeV -> TOT
	 double chargeCalib;     //MeV pr. eh-pair
	 int    chargeThreshold; //Threshold of preamp in #eh-pairs
	 double totThreshold;    // Same as above, but in TOT units
	 //   (calculated from above)
	 */
	double chargeThreshold; //Threshold of preamp in #eh-pairs
	double tot_a0, tot_a1, tot_a2;
	double tot_qMin, tot_qMax;
	double si_w;

	// Electrode parameters
	int nElec; //Number of readout electrodes pr. pixel
	double R_readout;
	double R_bias;

	//Efficiency inside bias and readout holes
	double eff_readout;
	double eff_bias;

	//Squared distance between elecPos and (xmod,ymod)
	inline double elecDist2(double elecPosX, double elecPosY, double xmod,
			double ymod) {
		return (elecPosX - xmod) * (elecPosX - xmod)
				+ (elecPosY - ymod) * (elecPosY - ymod);
	}

	//TOT response function
	inline int tot_response(double Q) {
		return (int) (tot_a0 * (tot_a1 + Q) / (tot_a2 + Q));
	}

public:
	Full3D_HP(const Event &event);
	virtual void init(TbConfig &config, const Event &event);
	virtual void digitize(Event& event, const TbConfig& config,
			std::vector<PllHit*>& hits);
	virtual void finalize(const TbConfig &config);
	~Full3D_HP();
};

#endif //FULL3D_HP_H
