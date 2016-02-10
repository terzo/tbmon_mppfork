#include "simDutEdep.h"

void simDutEdep::init(TbConfig &config) {
	if (not config.isSimulation) {
		//Todo: Rather check that edepBuilder is present.
		//      (ideally, edepBuilder etc. should only run if
		//       an analysis requiring it is enabled)
		cout << "[ simDutEdep ]; SimDutEdep not usable on real data" << endl;
		return;
	}

	h_sumEdep = new TH1D("", "Deposited energy in DUT", 400, 0, 0.30);
	h_sumEdep->SetXTitle("Deposited energy [MeV]");

	//config.fNeedSimEdep = true;
}

void simDutEdep::event(const TbConfig& config, const Event& event) {
	if (not config.isSimulation)
		return;
	if (event.fTrack != event::kGood || event.fSimSync != event::kGood)
		return;

	double sumEdep = 0.0;
	int incr = 0;
	for (vector<simPixelEdep>::iterator edep = event.simData->edeps.begin();
			edep != event.simData->edeps.end(); edep++) {

		sumEdep += (*edep).edep;
		incr++;
	}
	if (sumEdep > 0.0) {
		h_sumEdep->Fill(sumEdep);
	}

	/*
	 if (config.currentEntry % 100 == 0 && event.dut->iden == 160)
	 cout << config.currentEntry << " " << incr << " " << sumEdep<< endl;
	 */
}

void simDutEdep::finalize(const TbConfig& config) {
	if (not config.isSimulation)
		return;

	config.drawToFile(this->name, "sumEdep", h_sumEdep);
	config.saveToFile(this->name, "sumEdep", h_sumEdep);
}
