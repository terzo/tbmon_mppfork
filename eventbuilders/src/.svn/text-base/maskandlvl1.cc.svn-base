#include "maskandlvl1.h"

void MaskAndLvl1::init(TbConfig& config) {
	string lvl1string = config.cmdLineExtras_argGetter(
			"B_maskandlvl1_lvl1values", string(""), "Lvl1 values");
	if (lvl1string == "")
		lvl1cut = false;
	else {
		lvl1cut = true;
		string ints;
		stringstream lvl1stream(lvl1string);
		while (getline(lvl1stream, ints, ' ')) {
			lvl1.push_back(atoi(ints.c_str()));
		}
		//TBALOG (kDEBUG) << "Reading in lvl1 cut values: ";
		//for(vector<int>::const_iterator i = lvl1.begin(); i != lvl1.end(); i++) {TBALOG (kDEBUG) << (*i) << " ";}
	}
}

void MaskAndLvl1::buildEvent(Event &event, map<int, Event>&, TbConfig &config) {
	for (vector<PllHit*>::iterator it = event.rawHits.begin();
			it != event.rawHits.end(); it++) {
		int col((*it)->col), row((*it)->row);
		bool masked = false;
		bool cut = lvl1cut ? true : false;
		for (vector<pair<int, int> >::iterator mm = event.dut->masks.begin();
				mm != event.dut->masks.end(); mm++) {

			if ((*mm).first == row and (*mm).second == col) {
				masked = true;
				break;
			}
		}
		if (lvl1cut) {
			int lv1((*it)->lv1);
			for (vector<int>::const_iterator i = lvl1.begin(); i != lvl1.end();
					i++) {
				if (lv1 == (*i)) {
					cut = false;
					break;
				}
			}
		}
		if ((not masked) && (not cut)) {
			event.hits.push_back((*it));
		}
	}
}

