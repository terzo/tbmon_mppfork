#ifndef CHI2BUILDER_H
#define CHI2BUILDER_H

// standard header files
#include <vector>
#include <sstream>
#include <stdlib.h>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief Applies Chi2 cut by settings flags
 *
 */
class Chi2Builder: public EventBuilder {
private:
	double cutVal;
	double chi2default;
public:
	//This builder gets its default value from constructor argument
	Chi2Builder(double chi2) :
			chi2default(chi2) {
		name = "Chi2Builder";
	}
	;
	virtual void init(TbConfig& config);
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
};

#endif //CHI2BUILDER_H
