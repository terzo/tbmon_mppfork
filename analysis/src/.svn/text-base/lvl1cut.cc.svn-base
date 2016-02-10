#include "lvl1cut.h"

void Lvl1Cut::init(TbConfig &config) {

	h_lvl1HistAny = new TH1I("", ";All Hits Lvl1", 16, 0, 16);
	h_lvl1HistMatched = new TH1I("", ";Matched Hits Lvl1", 16, 0, 16);

}

void Lvl1Cut::event(const TbConfig &config, const Event &event) {

	for (vector<PllHit*>::const_iterator i = event.hits.begin();
			i != event.hits.end(); i++) {
		h_lvl1HistAny->Fill((*i)->lv1);
	}

	int matched = cluster::getMatched(event.clusters, event);
	if (matched == -1)
		return;

	for (vector<PllHit*>::const_iterator i = event.clusters.at(matched).begin();
			i != event.clusters.at(matched).end(); i++) {
		h_lvl1HistMatched->Fill((*i)->lv1);
	}

}

void Lvl1Cut::finalize(const TbConfig &config) {

	config.drawAndSave(name, "lvl1distany", h_lvl1HistAny);
	config.drawAndSave(name, "lvl1distmatched", h_lvl1HistMatched);

}
