#include "dutsync.h"

void DutSync::readFile(const char* fileName) {
	if (fileName == NULL) {
		return;
	}
	fstream stream;
	stream.open(fileName);
	if (not stream.is_open()) {
		cout << "Unable to open file: " << fileName << endl;
		exit(-1);
	}
	while (not stream.eof()) {
		int run(-1);
		int entries(0);
		stream >> run;
		stream >> entries;
		vector<bool> syncVector;
		if (entries == 1) {
			syncVector.push_back(true);
			syncMap[run] = syncVector;
			continue;
		}
		for (int ii = 0; ii < entries; ii++) {
			int accept;
			stream >> accept;
			if (accept == 1) {
				syncVector.push_back(true);
			} else if (accept == 0) {
				syncVector.push_back(false);
			} else {
				cout << "[ DUTSync ]:Error in stream, expected 0 or 1" << endl;
				exit(-1);
			}
		}
		syncMap[run] = syncVector;
	}
	int dummy;
	stream >> dummy;
}

void DutSync::initRun(TbConfig& config) {
	acceptSync = syncMap[config.currentRun];
	if (acceptSync.size() == 1 and acceptSync.at(0) == true) {
		checkCurrent = false;
	} else {
		checkCurrent = true;
	}
}

void DutSync::buildEvent(Event &event, map<int, Event> &events,
		TbConfig& config) {
	if (config.logLevel >= kDEBUG3)
		cout << "[ DUTSync ]: in buildEvent." << endl;

	if (checkCurrent == false) {
		event.fDutSync = event::kGood;
	} else {
		if (acceptSync.at(config.currentEntry) == 1) {
			event.fDutSync = event::kGood;
		} else {
			if (config.logLevel >= kDEBUG) {
				cout << "[ DUTSync ]; Event " << config.currentEntry
						<< " appears to be out of sync, marking track as bad."
						<< endl;
			}
			event.fDutSync = event::kBad;
			event.fTrack = event::kBad;
		}
	}
}
