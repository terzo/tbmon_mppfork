#include "batcheckdutsync.h"

void BatCheckDUTSync::falseSinceMatch(const TbConfig& config) {
	for (int ii = lastSurePoint; ii < syncVec.size(); ii++) {
		syncVec.at(ii) = false;
	}
	//We are not sure if these are bad or good.
	lastSurePoint = config.currentEntry;
}

void BatCheckDUTSync::init(TbConfig &config) {
	char* streamName = config.getOutStreamName(this->name, "dutsync");
	TBALOG(kINFO) << "Output to " << streamName << endl;
	stream.open(streamName);
	delete[] streamName;
}

void BatCheckDUTSync::initRun(const TbConfig &config) {
	syncVec.clear();
	syncVec.resize(config.tree->GetEntries());
	for (int ii = 0; ii < syncVec.size(); ii++) {
		syncVec.at(ii) = false;
	}
	assumeSync = true;
	perfect = true;
	lastSurePoint = 0;
	deltaTrig = 0;
}

void BatCheckDUTSync::event(const TbConfig &config, const Event &event) {
	//If there are no DUT hits there is really nothing to do
	if (event.track == NULL) {
		TBALOG(kERROR) <<
				"This analysis is intended for BAT data, but cannot find Track object."
				<< endl;
		return;
	}
	if (event.hits.size() < 1) {
		syncVec.at(config.currentEntry) = assumeSync;
		return;
	}
	//We have hits, lets check the delta trig
	int tmpDelta = event.track->trig - event.hits.at(0)->trig;
	if (tmpDelta != deltaTrig) {
		deltaTrig = tmpDelta;
		TBALOG(kDEBUG) << "Delta trig changed(" << deltaTrig
				<< "), sync is assumed lost at: " << config.currentEntry
				<< " Whiping all events back to last assumed known state: "
				<< lastSurePoint << endl;
		//We have had a sync change. We can no longer assume to be in sync
		assumeSync = false;
		perfect = false;
		//We are note sure of the state since the last sure point
		falseSinceMatch(config);
	}
	//We cannot be sure about a match if we do not trust the track
	if (event.fTrackAngle != event::kGood or event.fTrackChi2) {
		syncVec.at(config.currentEntry) = assumeSync;
		return;
	}
	//do we have a match?
	if (cluster::isMatch(event.hits, event)) {
		//Assume we are in sync
		lastSurePoint = config.currentEntry;
		if (assumeSync == false) {
			TBALOG(kDEBUG) << "Found match, assume we are back in sync: "
					<< config.currentEntry << endl;
		}
		assumeSync = true;
	}
	syncVec.at(config.currentEntry) = assumeSync;
}

void BatCheckDUTSync::finalizeRun(const TbConfig &config) {
	falseSinceMatch(config);
	stream << config.currentRun << " ";
	if (perfect) {
		stream << "1 ";
		return;
	}
	stream << syncVec.size() << " ";
	for (int ii = 0; ii < syncVec.size(); ii++) {
		if (syncVec.at(ii)) {
			stream << "1 ";
		} else {
			stream << "0 ";
		}
	}
}

void BatCheckDUTSync::finalize(const TbConfig &config) {
	stream.close();
}
