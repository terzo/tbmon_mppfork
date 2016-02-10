#include "chi2builder.h"

void Chi2Builder::init(TbConfig& config) {
	cutVal = config.cmdLineExtras_argGetter("G_chi2cut", chi2default,
			"Mark tracks with chi2 > this value as bad");
	if (cutVal <= 0.0) {
		cout << "[ Chi2Builder ]; Error: chi2 must be cut at value > 0" << endl;
	}
}

void Chi2Builder::buildEvent(Event &event, map<int, Event>&, TbConfig &) {
	//if(event.track->ct > cutVal){ // changed by BJD
	if (event.chi2 > cutVal) {
		event.fTrackChi2 = event::kBad;
		event.fTrack = event::kBad;
	} else {
		event.fTrackChi2 = event::kGood;
	}
}
