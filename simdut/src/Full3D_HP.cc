#include "Full3D_HP.h"

using namespace std;

Full3D_HP::Full3D_HP(const Event &event) {
	modelName = "Full3D_HP";

	//Make a memory-contiguous matrix
	qMatrix = new double*[event.dut->getNrows()];
	double* qMatrixRaw = new double[event.dut->getNcols()
			* event.dut->getNrows()];
	for (int row = 0; row < event.dut->getNrows(); row++) {
		qMatrix[row] = &qMatrixRaw[row * event.dut->getNcols()];
	}

	readoutE_xpos = NULL;
	biasE_xpos = NULL;
}

Full3D_HP::~Full3D_HP() {
	//delete &qMatrix; //This *should* get the whole qMatrixRaw // ?!?
	delete qMatrix[0]; //Pointer to the first element in qMatrixRaw. This should get the whole qMatrixRaw
	delete qMatrix; //Kill the "spine"

	if (readoutE_xpos != NULL)
		delete readoutE_xpos;
	if (biasE_xpos != NULL)
		delete biasE_xpos;
}

void Full3D_HP::init(TbConfig &config, const Event &event) {
	//Default values
	/* OLD model
	 totCalib = 60/0.069;
	 chargeCalib = 3.62e-6;
	 chargeThreshold = 3200;
	 totThreshold = chargeThreshold*chargeCalib*totCalib;
	 */
	chargeThreshold = 3200;

	// TOT fit from STA-3E
	tot_a0 = 731.521;
	tot_a1 = -1134.26;
	tot_a2 = 230466;
	tot_qMax = 45000;
	tot_qMin = 4000;

	si_w = 3.62e-6;  // MeV / EHP

	nElec = 3;
	R_readout = 7.0;
	R_bias = 7.0;

	//Calculate readout and bias electrode positions
	readoutE_xpos = new double[nElec];
	for (int i = 0; i < nElec; i++) {
		readoutE_xpos[i] = event.dut->getPitchX()
				* (i / (double) nElec - 0.5 + 1 / (double) (2.0 * nElec));
	}
	biasE_xpos = new double[nElec + 1];
	for (int i = 0; i < nElec + 1; i++) {
		biasE_xpos[i] = event.dut->getPitchX() * (i / (double) nElec - 0.5);
	}

	//Read from config.cmdLineExtras
	R_readout = config.cmdLineExtras_argGetter("SD_Full3D_HP_R_readout",
			R_readout, "Radius of readout electrodes [um]");
	R_bias = config.cmdLineExtras_argGetter("SD_Full3D_HP_R_bias", R_bias,
			"Radius of bias electrodes    [um]");

	eff_readout = config.cmdLineExtras_argGetter("SD_Full3D_HP_eff_readout",
			0.0, "Collection efficiency inside readout electrode");
	eff_bias = config.cmdLineExtras_argGetter("SD_Full3D_HP_eff_bias", 0.0,
			"Collection efficiency inside bias electrode");
}

void Full3D_HP::digitize(Event& event, const TbConfig& config,
		std::vector<PllHit*>& hits) {

	//Zero qMatrix
	for (int row = 0; row < event.dut->getNrows(); row++) {
		for (int col = 0; col < event.dut->getNcols(); col++) {
			qMatrix[row][col] = 0.0;
		}
	}

	double si_w_inv = 1.0 / si_w;

	//Fill qMatrix
	for (vector<simPixelEdep>::iterator edep = event.simData->edeps.begin();
			edep != event.simData->edeps.end(); edep++) {

		//Coordinates with origo in center of pixel (0,0) (col, row)
		double xp = (*edep).posLocal.data[0];
		double yp = (*edep).posLocal.data[1];
		//Row and column of current pixel
		int col = (int) ((xp + event.dut->getPitchX() / 2.0)
				/ event.dut->getPitchX());
		int row = (int) ((yp + event.dut->getPitchY() / 2.0)
				/ event.dut->getPitchY());
		//Coordinates local to the pixel we are in, origo at pixel center
		double xmod = xp - col * event.dut->getPitchX();
		double ymod = yp - row * event.dut->getPitchY();

		//check if we are inside hole
		bool insideHole = false;
		bool insideBias = false;
		for (int i = 0; i < nElec; i++) { //Readout
			if (elecDist2(readoutE_xpos[i], 0.0, xmod, ymod)
					< R_readout * R_readout) {
				insideHole = true;
				break;
			}
		}
		if (!insideHole) {
			for (int i = 0; i < nElec + 1; i++) { //Bias
				if (elecDist2(biasE_xpos[i], -event.dut->getPitchY() / 2, xmod,
						ymod) < R_bias * R_bias
						|| elecDist2(biasE_xpos[i], +event.dut->getPitchY() / 2,
								xmod, ymod) < R_bias * R_bias) {
					insideHole = true;
					insideBias = true;
					break;
				}
			}
		}
		//if (insideHole) break; //No edep

		//Possible to get edep at exact edge of sensor
		if (row == event.dut->getNrows())
			row--;
		if (col == event.dut->getNcols())
			col--;
		if (!insideHole) {
			//Normal bulk
			qMatrix[row][col] += si_w_inv * (*edep).edep;
		} else if (!insideBias) {
			//Inside readout
			qMatrix[row][col] += si_w_inv * (*edep).edep * eff_readout;
		} else {
			//Inside bias
			qMatrix[row][col] += si_w_inv * (*edep).edep * eff_bias;
		}
	}

	//Generate digits
	for (int row = 0; row < event.dut->getNrows(); row++) {
		for (int col = 0; col < event.dut->getNcols(); col++) {
			if (qMatrix[row][col] > chargeThreshold) {
				PllHit* hit = new PllHit();
				hit->trig = event.track->trig;
				hit->iden = event.dut->getDUTid();
				hit->chp = 0;
				hit->row = row;
				hit->col = col;
				hit->tot = tot_response(qMatrix[row][col]);
				hit->lv1 = 0;
				hits.push_back(hit);

			}
		}
	}
}

void Full3D_HP::finalize(const TbConfig &config) {
	int N = 50;
	double* Q = new double[N];
	double Qh = (tot_qMax - chargeThreshold) / ((double) (N - 1));
	double* ToT = new double[N];
	for (int i = 0; i < N; i++) {
		Q[i] = chargeThreshold + Qh * i;
		ToT[i] = tot_response(Q[i]);
	}

	TGraph* g_resp = new TGraph(N, Q, ToT);
	g_resp->SetName("Full3D_HP::resp");
	g_resp->SetTitle("Full3D_HP::resp");

	config.drawToFile(this->modelName.c_str(), "resp", "AL", g_resp);
	config.saveToFile(this->modelName.c_str(), "resp", g_resp);

}

