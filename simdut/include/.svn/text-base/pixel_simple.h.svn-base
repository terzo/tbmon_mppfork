#ifndef PIXEL_SIMPLE_H
#define PIXEL_SIMPLE_H

// tbmon header files
#include "simdut.h"
#include "event.h"
#include "simDutEdep.h"

/*! \brief This is a very simple implementation of a simulation pixel model.
 *
 * This toy model is 100% efficient everywhere,
 * with TOT response proportional to edep,
 * and always LVL1 = 1.
 */
class pixel_simple: public Simdut {
private:
	double** totMatrix; //Non-sparse storage of 'digits' TOT values

	double totCalib;        //Conversion factor from MeV -> TOT
	double chargeCalib;     //MeV pr. eh-pair
	int chargeThreshold; //Threhsold of preamp in #eh-pairs
	double totThreshold;    // Same as above, but in TOT units
							//   (calculated from above)
public:
	pixel_simple(const Event &event);
	virtual void init(TbConfig &config);
	virtual void digitize(Event& event, const TbConfig& config,
			std::vector<PllHit*>& hits);

};

#endif //PIXEL_SIMPLE_H
