#include "Full3D_Vadim.h"

using namespace std;

Full3D_Vadim::Full3D_Vadim(const Event &event) {
	modelName = "Full3D_Vadim";

	qMatrix = new double*[event.dut->getNrows()];
	double* qMatrixRaw = new double[event.dut->getNcols()
			* event.dut->getNrows()];
	for (int row = 0; row < event.dut->getNrows(); row++) {
		qMatrix[row] = &qMatrixRaw[row * event.dut->getNcols()];
	}

	readoutE_xpos = NULL;
	biasE_xpos = NULL;
}

Full3D_Vadim::~Full3D_Vadim() {
	//delete &totMatrix; //This *should* get totMatrixRaw ??
	delete qMatrix[0]; //Pointer to the first element in qMatrixRaw. This should get the whole qMatrixRaw
	delete qMatrix; //Kill the "spine"

	if (readoutE_xpos != NULL)
		delete readoutE_xpos;
	if (biasE_xpos != NULL)
		delete biasE_xpos;
}

void Full3D_Vadim::init(TbConfig &config, const Event &event) {

	//Set defaults and get arguments
	/*
	 //Old response model
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
	cout << "Getting configs..." << endl << flush;
	R_readout = config.cmdLineExtras_argGetter("SD_Full3D_Vadim_R_readout", 5.0,
			"Radius of readout electrodes                [um]");
	R_bias = config.cmdLineExtras_argGetter("SD_Full3D_Vadim_R_bias", 5.0,
			"Radius of bias electrodes                   [um]");
	s_readout = config.cmdLineExtras_argGetter("SD_Full3D_Vadim_S_readout", 0.5,
			"Efficiency fall-off distance around readout [um]");
	s_bias = config.cmdLineExtras_argGetter("SD_Full3D_Vadim_S_bias", 0.5,
			"Efficiency fall-off distance around bias    [um]");

	cout << "Done." << endl << flush;

	//Calculate readout and bias electrode positions (2x2 pixel cell)
	if (nElec % 2 == 0) {
		//Even number of readout electrodes pr. pixel

		double delta = 1 / (double) nElec; //Inter-electrode distance

		N_readout = 2 * nElec;
		readoutE_xpos = new double[N_readout];
		for (int i = 0; i < N_readout; i++) {
			readoutE_xpos[i] = (-1.0 + delta * (i + 0.5))
					* event.dut->getPitchX();
		}

		N_bias = 2 * nElec + 1;
		biasE_xpos = new double[N_bias];
		for (int i = 0; i < N_bias; i++) {
			biasE_xpos[i] = (-1.0 + delta * i) * event.dut->getPitchX();
		}
	} else {
		//Odd number of electrodes

		double delta = 1 / (double) nElec; //Inter-electrode distance

		N_readout = nElec + 2 * (int) ceil(nElec / 2);
		readoutE_xpos = new double[N_readout];
		for (int i = 0; i < N_readout; i++) {
			readoutE_xpos[i] = (-1.0 + delta * i) * event.dut->getPitchX();
		}

		N_bias = N_readout - 1;
		biasE_xpos = new double[N_bias];
		for (int i = 0; i < N_bias; i++) {
			biasE_xpos[i] = (-1.0 + delta * (i + 0.5)) * event.dut->getPitchX();
		}
	}
}

void Full3D_Vadim::digitize(Event& event, const TbConfig& config,
		std::vector<PllHit*>& hits) {

	double si_w_inv = 1.0 / si_w;

	//Zero qMatrix
	for (int row = 0; row < event.dut->getNrows(); row++) {
		for (int col = 0; col < event.dut->getNcols(); col++) {
			qMatrix[row][col] = 0.0;
		}
	}

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

		//Electrode efficiency
		double eff = 1.0;
		double r = 0.0;

		//Readout
		for (int i = 0; i < N_readout; i++) {
			r = sqrt(elecDist2(readoutE_xpos[i], 0.0, xmod, ymod));
			eff /= 1 + exp((R_readout - r) / s_readout);

			r = sqrt(
					elecDist2(readoutE_xpos[i], event.dut->getPitchY(), xmod,
							ymod));
			eff /= 1 + exp((R_readout - r) / s_readout);

			r = sqrt(
					elecDist2(readoutE_xpos[i], -event.dut->getPitchY(), xmod,
							ymod));
			eff /= 1 + exp((R_readout - r) / s_readout);
		}
		//Bias
		for (int i = 0; i < N_bias; i++) {
			r = sqrt(
					elecDist2(biasE_xpos[i], 0.5 * event.dut->getPitchY(), xmod,
							ymod));
			eff /= 1 + exp((R_bias - r) / s_bias);

			r = sqrt(
					elecDist2(biasE_xpos[i], -0.5 * event.dut->getPitchY(),
							xmod, ymod));
			eff /= 1 + exp((R_bias - r) / s_bias);
		}

		//Possible to get edep at exact edge of sensor: Drop this edep
		if (row == event.dut->getNrows() || col == event.dut->getNcols()) {
			if (config.logLevel == kDEBUG3) {
				cout
						<< "[ Full3D_Vadim ]::digitize(); Got chargeDep at sensor edge, dropping"
						<< endl;
			}
			continue;
		}

		//cout << "(" << xmod << "," << ymod "), eff=" << eff << endl;

		//totMatrix[row][col] += totCalib*(*edep).edep*eff;
		qMatrix[row][col] += si_w_inv * (*edep).edep * eff;
	}

	//Generate digits from totMatrix content
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

double Full3D_Vadim::elecDist2(double elecPosX, double elecPosY, double xmod,
		double ymod) {
	return (elecPosX - xmod) * (elecPosX - xmod)
			+ (elecPosY - ymod) * (elecPosY - ymod);
}

void Full3D_Vadim::finalize(const TbConfig& config) {
	//Make nice plot of electrode efficiencies with current parameters
	double r_max = R_bias > R_readout ? R_bias * 3 : R_readout * 3;
	int N = 100;
	double r_delta = r_max / (double) N;
	double* r = new double[N];
	double* eff_readout = new double[N];
	double* eff_bias = new double[N];
	for (int i = 0; i < N; i++) {
		r[i] = i * r_delta;
		eff_readout[i] = 1.0 / (1.0 + exp((R_readout - r[i]) / s_readout));
		eff_bias[i] = 1.0 / (1.0 + exp((R_bias - r[i]) / s_bias));
	}

	TGraph* g_eff_readout = new TGraph(N, r, eff_readout);
	g_eff_readout->SetName("Full3D_vadim::eff_readout");
	g_eff_readout->SetTitle("Full3D_vadim::eff_readout");
	TGraph* g_eff_bias = new TGraph(N, r, eff_bias);
	g_eff_bias->SetName("Full3D_vadim::eff_bias");
	g_eff_bias->SetTitle("Full3D_vadim::eff_bias");

	config.drawToFile(this->modelName.c_str(), "eff_readout", "AL",
			g_eff_readout);
	config.saveToFile(this->modelName.c_str(), "eff_readout", g_eff_readout);
	config.drawToFile(this->modelName.c_str(), "eff_bias", "AL", g_eff_bias);
	config.saveToFile(this->modelName.c_str(), "eff_bias", g_eff_bias);
}
