#include "anglecuts.h"

AngleCuts::AngleCuts(const char* calibFile, double sigmas) {
	name = "AngleCuts";

	this->sigmas = sigmas;
	if (calibFile == NULL) {
		return;
	}
	fstream stream;
	stream.open(calibFile);
	if (not stream.is_open()) {
		cerr << "Unable to open file: " << calibFile << endl;
		exit(-1);
	}
	//Dump calib filt to map
	while (not stream.eof()) {
		AngleStats stats;
		stream >> stats.run;
		stream >> stats.meanX;
		;
		stream >> stats.sigmaX;
		stream >> stats.meanY;
		stream >> stats.sigmaY;
		angleStats[stats.run] = stats;
	}
}
void AngleCuts::initRun(TbConfig &config) {
	//Get angle stats for current run
	if (angleStats.find(config.currentRun) == angleStats.end()) {
		std::cerr << "Can not find angle calibration for run"
				<< config.currentRun << ", will quit." << endl
				<< "Remove AngleCuts eventbuilder and rerun calibrations."
				<< endl;
		exit(-1);
	}
	this->currentCuts = angleStats[config.currentRun];
}

void AngleCuts::buildEvent(Event &event, map<int, Event> &, TbConfig &) {
	//Check angles in plane with iden 0
	double angleX(99.0), angleY(99.0);
	event.fTrackAngle = event::kGood;
	for (int par = 0; par < event.track->nTrackParams; par++) {
		TrackParams* params = (TrackParams*) event.track->trackParams[par];
		if (params->iden != 0) {
			continue;
		}
		angleX = params->params[2];
		angleY = params->params[3];
		break;
	}
	event.fTrackAngle = event::kGood;
	//If the angles are bad, so is the track
	if (fabs(angleX - currentCuts.meanX) > (sigmas * currentCuts.sigmaX)) {
		event.fTrackAngle = event::kBad;
		event.fTrack = event::kBad;
	}
	if (fabs(angleY - currentCuts.meanY) > (sigmas * currentCuts.sigmaY)) {
		event.fTrackAngle = event::kBad;
		event.fTrack = event::kBad;
	}
}
