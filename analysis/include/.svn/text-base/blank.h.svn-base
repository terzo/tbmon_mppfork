#ifndef BLANK_H
#define BLANK_H

// tbmon header files
#include "event.h"
#include "tbanalysis.h"

/*! \brief Empty skeleton analysis
 *
 */
class Blank: public TbAnalysis {

public:
	virtual void init(TbConfig &config);
	virtual void initRun(const TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalizeRun(const TbConfig &config);
	virtual void finalize(const TbConfig &config);
};

#endif //BLANK_H
