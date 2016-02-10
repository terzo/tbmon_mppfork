#ifndef BATCHECKDUTSYNCH_H
#define BATCHECKDUTSYNCH_H

// standard header files
#include <vector>
#include <fstream>
#include <iostream>

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief ONLY FOR BAT DATA (no longer maintained - will be probably removed at some point): checks synchronisation between DUTs and telescope
 *
 */
class BatCheckDUTSync: public TbAnalysis {
private:
	vector<bool> syncVec;
	ofstream stream;
	int lastSurePoint;
	int deltaTrig;
	bool assumeSync;
	bool perfect;
	void falseSinceMatch(const TbConfig &config);
	void dumpVectorToStream();
public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
	virtual void initRun(const TbConfig &config);
	virtual void finalizeRun(const TbConfig &config);
};

#endif //BATCHECKDUTSYNCH_H
