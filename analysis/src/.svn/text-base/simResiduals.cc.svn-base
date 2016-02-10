#include "simResiduals.h"

void simResiduals::init(TbConfig &config) {
	if (not config.isSimulation) {
		//Todo: Rather check that truthBuilder is present.
		//      (ideally, truthBuilder etc. should only run if
		//       an analysis requiring it is enabled)
		cout << "[ simResiduals ]; SimResiduals not usable on real data"
				<< endl;
		return;
	}
	h_resX = new TH1D("", "Track vs. truth residuals (X)", 400, -50, 50);
	h_resX->SetXTitle("Track - truth  [#mu m]");
	h_resY = new TH1D("", "Track vs. truth residuals (Y)", 400, -50, 50);
	h_resY->SetXTitle("Track - truth  [#mu m]");

	h_XvsX = new TH2D("", "TruthX vs TrackX correlation", 100, 0, 10000, 100, 0,
			10000);
	h_XvsY = new TH2D("", "TruthX vs TrackY correlation", 100, 0, 10000, 100, 0,
			10000);
	h_YvsX = new TH2D("", "TruthY vs TrackX correlation", 100, 0, 10000, 100, 0,
			10000);
	h_YvsY = new TH2D("", "TruthY vs TrackY correlation", 100, 0, 10000, 100, 0,
			10000);

	//config.fNeedSimTruth = true;
}

void simResiduals::event(const TbConfig& config, const Event& event) {
	if (not config.isSimulation)
		return;
	if (event.fTrack != event::kGood || event.fSimSync != event::kGood)
		return;
	if (event.simData->truthHits.size() != 1)
		return; //Drop events with multiple or zero truthHits

	h_resX->Fill(event.trackX - event.simData->truthHits[0]->posLocal.data[0]);
	h_resY->Fill(event.trackY - event.simData->truthHits[0]->posLocal.data[1]);
	/*
	 cout << event.trackX << "\t" << event.simData->truthHits[0]->posLocal.data[0] << endl;
	 cout << event.trackY << "\t" << event.simData->truthHits[0]->posLocal.data[1] << endl;
	 cout << endl;
	 */

	h_XvsX->Fill(event.simData->truthHits[0]->posLocal.data[0], event.trackX);
	h_XvsY->Fill(event.simData->truthHits[0]->posLocal.data[0], event.trackY);
	h_YvsX->Fill(event.simData->truthHits[0]->posLocal.data[1], event.trackX);
	h_YvsY->Fill(event.simData->truthHits[0]->posLocal.data[1], event.trackY);

}

void simResiduals::finalize(const TbConfig& config) {
	if (not config.isSimulation)
		return;

	config.drawToFile(this->name, "resX", h_resX);
	config.drawToFile(this->name, "resY", h_resY);
	config.saveToFile(this->name, "resX", h_resX);
	config.saveToFile(this->name, "resY", h_resY);

	config.drawToFile(this->name, "XvsX", h_XvsX);
	config.drawToFile(this->name, "XvsY", h_XvsY);
	config.drawToFile(this->name, "YvsX", h_YvsX);
	config.drawToFile(this->name, "YvsY", h_YvsY);
}
