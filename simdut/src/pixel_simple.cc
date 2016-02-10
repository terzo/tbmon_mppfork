#include "pixel_simple.h"

using namespace std;

pixel_simple::pixel_simple(const Event &event) {
	modelName = "pixel_simple";

	//Initialize storage
	totMatrix = new double*[event.dut->getNrows()];
	double* totMatrixRaw = new double[event.dut->getNcols()
			* event.dut->getNrows()];
	for (int row = 0; row < event.dut->getNrows(); row++) {
		totMatrix[row] = &totMatrixRaw[row * event.dut->getNcols()];
	}

}

void pixel_simple::init(TbConfig &config) {
	//Default values
	totCalib = 60 / 0.069;
	chargeCalib = 3.62e-6;
	chargeThreshold = 3200;
	totThreshold = chargeThreshold * chargeCalib * totCalib;

	//Read from config.cmdLineExtras
	//TODO
}

void pixel_simple::digitize(Event& event, const TbConfig& config,
		std::vector<PllHit*>& hits) {

	//Zero totMatrix
	for (int row = 0; row < event.dut->getNrows(); row++) {
		for (int col = 0; col < event.dut->getNcols(); col++) {
			totMatrix[row][col] = 0.0;
		}
	}

	//Fill totMatrix
	for (vector<simPixelEdep>::iterator edep = event.simData->edeps.begin();
			edep != event.simData->edeps.end(); edep++) {

		int col = (int) (((*edep).posLocal.data[0]
				+ event.dut->getPitchX() / 2.0) / event.dut->getPitchX());
		int row = (int) (((*edep).posLocal.data[1]
				+ event.dut->getPitchY() / 2.0) / event.dut->getPitchY());

		//Possible to get edep at exact edge of sensor
		if (row == event.dut->getNrows())
			row--;
		if (col == event.dut->getNcols())
			col--;

		totMatrix[row][col] += totCalib * (*edep).edep;
	}

	//Generate digits
	for (int row = 0; row < event.dut->getNrows(); row++) {
		for (int col = 0; col < event.dut->getNcols(); col++) {
			if (totMatrix[row][col] > totThreshold) {
				PllHit* hit = new PllHit();
				hit->trig = event.track->trig;
				hit->iden = event.dut->getDUTid();
				hit->chp = 0;
				hit->row = row;
				hit->col = col;
				hit->tot = (int) totMatrix[row][col];
				hit->lv1 = 0;
				hits.push_back(hit);

			}
		}
	}
}
