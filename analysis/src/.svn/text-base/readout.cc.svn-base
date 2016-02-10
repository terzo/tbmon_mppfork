#include "readout.h"

void Readout::init(TbConfig &config) {

}

void Readout::event(const TbConfig &config, const Event &event) {
	const char* eventflagnames[] = { "kBad", "kGood", "kUnknown" };

	//cout << "fBase:" << eventflagnames[event.fBase] << " ";
	//cout << "fClusters:" << eventflagnames[event.fClusters] << " "; // Are clusters built?
	//cout << "fTrack: " << eventflagnames[event.fTrack] << " ";
	//cout << "fTrackChi2:" << eventflagnames[event.fTrackChi2] << " "; // Is track Chi2 accepted?
	//cout << "fTrackAngle:" << eventflagnames[event.fTrackAngle] << " ";
	//cout << "fHits:" << eventflagnames[event.fHits] << " ";
	//cout << "fTrackCentralRegion:" << eventflagnames[event.fTrackCentralRegion] << " "; 
	//cout << "fTrackMaskedRegion:" << eventflagnames[event.fTrackMaskedRegion] << " "; 
	//cout << "fTrackRegion:" << eventflagnames[event.fTrackRegion] << " "; 
	//cout << "fTrackMatchNeighbour:" << eventflagnames[event.fTrackMatchNeighbour] << " ";

	//cout << endl << "Cluster size: " << (int) event.clusters.size() << endl;
	//event::eventflag_t fEtaCorrections;
	//event::eventflag_t fDutSync;
	//event::eventflag_t fSimSync; //Is simulation and reconstructed data well synced up?
	//event::eventflag_t fEtaCut;  //BAT eta cuts
	if (event.fTrackCentralRegion != event::kGood)
		return;
	if ((event.trackX < 0.0) || (event.trackY < 0.0)) {
		cout << "Negative" << endl;
		cout
				<< fmod(event.trackX + 0.5 * event.dut->getPitchX(),
						2 * event.dut->getPitchX()) << "   "
				<< fmod(event.trackY + 0.5 * event.dut->getPitchY(),
						2 * event.dut->getPitchY()) << endl;

		cout << endl;
	}
}

void Readout::finalize(const TbConfig &config) {
}

