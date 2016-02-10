#include "batangledist.h"

void BatAngleDist::init(TbConfig &config) {
	TBALOG(kDEBUG) << "Initializing histos...";
	char* sname = config.getOutStreamName(name, "trackangle");
	stream.open(sname);
	h_angleXTmp = new TH1D("", "", 1000, -0.01, 0.01);
	h_angleYTmp = new TH1D("", "", 1000, -0.01, 0.01);
	h_acceptedAnglesX = new TH1D("", "", 1000, -0.01, 0.01);
	h_rejectedAnglesX = new TH1D("", "", 1000, -0.01, 0.01);
	h_acceptedAnglesY = new TH1D("", "", 1000, -0.01, 0.01);
	h_rejectedAnglesY = new TH1D("", "", 1000, -0.01, 0.01);
	TBALOG(kDEBUG) << "Done" << endl;
}

void BatAngleDist::initRun(const TbConfig &config) {
	h_angleXTmp->Reset();
	h_angleYTmp->Reset();
}

void BatAngleDist::event(const TbConfig &config, const Event &event) {
	if (event.track == NULL) {
		TBALOG(kERROR)
				<< "This analysis is intended for BAT data, but cannot find Track object."
				<< endl;
		return;
	}
	TBALOG(kDEBUG3) << "In event" << endl;
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

void BatAngleDist::finalizeRun(const TbConfig &config) {
	stream << config.currentRun << " " << h_angleXTmp->GetMean() << " "
			<< h_angleXTmp->GetRMS() << " " << h_angleYTmp->GetMean() << " "
			<< h_angleYTmp->GetRMS() << " " << endl;
}

void BatAngleDist::finalize(const TbConfig &config) {
	stream.close();
	config.drawToFile(name, "acceptedX", h_acceptedAnglesX);
	config.drawToFile(name, "acceptedY", h_acceptedAnglesY);
	config.drawToFile(name, "rejectedX", h_rejectedAnglesX);
	config.drawToFile(name, "rejectedY", h_rejectedAnglesY);
}

