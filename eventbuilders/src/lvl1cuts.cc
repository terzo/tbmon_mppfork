#include "lvl1cuts.h"

void LVL1Cuts::buildEvent(Event &event, map<int, Event>&, TbConfig &config) {

	for (vector<PllHit*>::iterator it = event.rawHits.begin();
			it != event.rawHits.end(); it++) {
		const DUT* dut = config.getDut((*it)->iden);
		//Extract lv1 configuration from config
		int lv1Min = dut->lv1Min;
		int lv1Max = dut->lv1Max;
		int lv1((*it)->lv1);
//    TBALOG(kDEBUG3) << "Cuts from DUT calib: " << lv1Min << " < lv1 < " << lv1Max << endl;
		// filter hits with lvl1 out of range
		if (lv1Min <= lv1 || lv1 <= lv1Max)
			event.hits.push_back((*it));

	}
}

