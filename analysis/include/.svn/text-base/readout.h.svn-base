#ifndef READOUT_H
#define READOUT_H

#include "event.h"
#include "tbanalysis.h"

/*! \brief Seems sto be some test output of track coordinates.
 *
 */
class Readout: public TbAnalysis {

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //READOUT_H
