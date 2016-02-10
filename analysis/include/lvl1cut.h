#ifndef LVL1CUT_H
#define LVL1CUT_H

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Plots matched and un-matched hit lvl1 distributions.
 *
 */
class Lvl1Cut: public TbAnalysis {

private:
	TH1I* h_lvl1HistAny;
	TH1I* h_lvl1HistMatched;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);

};

#endif //LVL1CUT_H
