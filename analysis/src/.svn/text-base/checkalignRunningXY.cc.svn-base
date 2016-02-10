#include "checkalignRunningXY.h"
#include <cstring>
#include <string>

void CheckAlignRunningXY::init(TbConfig &config) {

}

void CheckAlignRunningXY::initRun(const TbConfig &config) {
	const DUT* dut = config.getDut(this->iden);
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();
	doCuts = true;
	_postShiftSubset = 1000;
	_firstPostShiftEvent = -1;
	_postShiftEventSkip = 50;

	_runningMeanX = new SensorInfoObjectCollection();
	_runningMeanY = new SensorInfoObjectCollection();

	_postShift = new PostShiftAlign();

	char* streamName = config.getOutStreamName(name, "translations");

	std::string originalStreamName = streamName;

	std::string firstStreamNamePart = originalStreamName.substr(0,originalStreamName.find(".")-1);

	std::stringstream runNumberString;
	runNumberString << config.currentRun;

	std::string newFileName = firstStreamNamePart + "-" + runNumberString.str() + ".root";

	_postShiftFile = new TFile(streamName,"RECREATE");

	_postShift->setDirectory(_postShiftFile);

	avX = new RunningAverage(dut->getDUTid(), _postShiftSubset);
	avY = new RunningAverage(dut->getDUTid(), _postShiftSubset);

	_runningMeanX->add(avX);
	_runningMeanY->add(avY);
}

void CheckAlignRunningXY::event(const TbConfig &config, const Event& event) {
	TBALOG(kDEBUG3) << "Started analysis class CheckAlignRunningXY event() command" << endl;
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		TBALOG (kDEBUG3) << "Event not passed all cuts" << endl;
		return;
	}
	// Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		TBALOG (kDEBUG3) << "Clusters not present in event" << endl;
		return;
	}
	// Are we in a good region of the chip?
	if (event.fTrackRegion != event::kGood) {
		TBALOG (kDEBUG3) << "Track not in good region" << endl;
		return;
	}

	int eventNum = config.currentEntry;
	int sensorId = config.getDut(this->iden)->getDUTid();

	// Look for matched clusters
	if (cluster::getMatched(event.clusters, event) != -1) {
		TBALOG(kDEBUG2) << "Found matching cluster!" << endl;
	} else return;

	float distX(0), distY(0);

	distX = event.trackX - cluster::getChargeWeightedCol(event.clusters[cluster::getMatched(event.clusters,event)]);
	distY = event.trackY - cluster::getChargeWeightedRow(event.clusters[cluster::getMatched(event.clusters,event)]);
	RunningAverage* avX = (RunningAverage*) (_runningMeanX->get(sensorId));
	RunningAverage* avY = (RunningAverage*) (_runningMeanY->get(sensorId));

	avX->add(distX);
	avY->add(distY);

	//if (eventNum % 10000 == 0) cout << "Queue fill level for sensor " << sensorId << ": " << avX->getFillRatio() << endl;
	double shiftX, shiftY;
	if (avX->isFilledUp() && avY->isFilledUp()) {
		if (!avX->isMeanBeginSet() && !avY->isMeanBeginSet()) {
			avX->setMeanBegin(eventNum / 2);
			avY->setMeanBegin(eventNum / 2);
			cout << "PostShift Queue filled for sensor " << sensorId
					<< " at event " << eventNum / 2 << endl;
			for (int i = 0; i < eventNum / 2; i++) {
				if ((eventNum / 2 - i) % _postShiftEventSkip == 0) {
					shiftX = avX->getMeanBegin(i);
					shiftY = avY->getMeanBegin(i);
					_postShift->add(sensorId, shiftX, shiftY, eventNum / 2 - i,
							_postShiftSubset);
				}
			}
		}
		shiftX = avX->getMean();
		shiftY = avY->getMean();
		int offset = avX->getMeanBegin(); // X and Y collection are identical for this
		// Is this is the collection for the first 1000 Hits the mean should be stored for hit 500
		if (eventNum % _postShiftEventSkip == 0)
			_postShift->add(sensorId, shiftX, shiftY, eventNum - offset,
					_postShiftSubset);
	}
}

void CheckAlignRunningXY::finalizeRun(const TbConfig &config) {
	_postShiftFile->Close();
	delete avX;
	delete avY;
	delete _postShift;
	delete _runningMeanX;
	delete _runningMeanY;
}

void CheckAlignRunningXY::finalize(const TbConfig &config) {

}
