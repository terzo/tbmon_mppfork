#include "batcorrelation.h"

void BatCorrelation::init(TbConfig &config, const Event &event) {
	= new TProfile("","",620,0,620);
	= new TProfile("","",620,0,620);
	= new TProfile("","",620,0,620);
	= new TProfile("","",620,0,620);
	= new TProfile("","",620,0,620);
	= new TProfile("","",620,0,620);
	char* sname = config.getOutStreamName(name, "batcorrelation");
	stream.open(sname);
}

void BatCorrelation::initRun(const TbConfig &config) {
	;
}

void BatCorrelation::event(const TbConfig &config, const Event &event) {
	if (event.track == NULL) {
		TBALOG(kERROR)
				<< "This analysis is intended for BAT data, but cannot find Track object."
				<< endl;
		return;
	}
	if (event.fTrackChi2 != event::kGood) {
		return;
	}
	double angleX(-11), angleY(-11);
	for (int ii = 0; ii < event.track->nTrackParams; ii++) {
		TrackParams* params = (TrackParams*) event.track->trackParams[ii];
		if (params->iden != 0) {
			continue;
		}
		angleX = params->params[2];
		angleY = params->params[3];
		break;
	}
	if (angleX == -11 or angleY == -11) {
		TBALOG(kINFO) << "Did not find track angles" << endl;
		return;
	}
	if (event.fTrackAngle == event::kGood) {
		h_acceptedAnglesX->Fill(angleX);
		h_acceptedAnglesY->Fill(angleY);
	}
	if (event.fTrackAngle == event::kBad) {
		h_rejectedAnglesX->Fill(angleX);
		h_rejectedAnglesY->Fill(angleY);
	}
	h_angleXTmp->Fill(angleX);
	h_angleYTmp->Fill(angleY);
}

void BatCorrelation::finalizeRun(const TbConfig &config) {
	stream << config.currentRun << " " << h_angleXTmp->GetMean() << " "
			<< h_angleXTmp->GetRMS() << " " << h_angleYTmp->GetMean() << " "
			<< h_angleYTmp->GetRMS() << " " << endl;
}

void BatCorrelation::finalize(const TbConfig &config, const Event &event) {
	stream.close();
	config.drawToFile(name, "acceptedX", h_acceptedAnglesX);
	config.drawToFile(name, "acceptedY", h_acceptedAnglesY);
	config.drawToFile(name, "rejectedX", h_rejectedAnglesX);
	config.drawToFile(name, "rejectedY", h_rejectedAnglesY);
}
