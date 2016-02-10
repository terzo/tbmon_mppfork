#include "battrack.h"

void BatTrack::initRun(TbConfig& config) {
	tracktree = (TTree*) config.infile->Get("Track");
	m_track = new Track();
	tracktree->SetBranchAddress("track", &m_track);
}

void BatTrack::finalizeRun(TbConfig& config) {
	if (not (m_track == NULL)) {
		delete m_track;
	}
	if (not (tracktree == NULL)) {
		delete tracktree;
	}
}
void BatTrack::initEvent(TbConfig &config, map<int, Event> &events) {
	tracktree->GetEntry(config.currentEntry);
}

void BatTrack::buildEvent(Event &event, map<int, Event> &, TbConfig & config) {
	bool found(false);
	event.track = m_track;
	for (int par = 0; par < event.track->nTrackParams; par++) {
		TrackParams* params = (TrackParams*) event.track->trackParams[par];
		if (params->iden != event.dut->getDUTid()) {
			continue;
		}
		found = true;
		event.setEvent(params->params[0] * 1000.0, params->params[1] * 1000.0,
				params->params[2], params->params[3], event.track->ct,
				event.track->ft);
		break;
	}
	event.fTrack = event::kGood;
	if (not found) {
		cout << "TbConfig::buildEvent: Track not found!" << endl;
		event.fBase = event::kBad;
	}
	//Adding pllHits to events
	for (int ii = 0; ii < event.track->nPllHit; ++ii) {
		PllHit* hit = (PllHit*) event.track->pllHit[ii];
		if (event.dut->getDUTid() == hit->iden) {
			event.rawHits.push_back(hit);
		}
	}
}
