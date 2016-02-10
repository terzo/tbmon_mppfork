#include "checkcentralregion.h"

void CheckCentralRegion::buildEvent(Event &event, map<int, Event> &,
		TbConfig & config) {

	//Check if track is inside the central region of the sensor

	double trackX, trackY, pitchX, pitchY;
	trackX = event.trackX - 0.5 * event.dut->getPitchX();
	trackY = event.trackY - 0.5 * event.dut->getPitchY();
	pitchX = event.dut->getPitchX();
	pitchY = event.dut->getPitchY();

	event.fTrackCentralRegion = event::kGood;

	if (trackX < pitchX * (event.dut->getSkipCols() - 1)) {
		event.fTrackCentralRegion = event::kBad;
	} else if (pitchX * (event.dut->getNcols() - (event.dut->getSkipCols() + 1))
			< trackX) {
		event.fTrackCentralRegion = event::kBad;
	} else if (trackY < pitchY * (event.dut->getSkipRows() - 1)) {
		event.fTrackCentralRegion = event::kBad;
	} else if (pitchY * (event.dut->getNrows() - (event.dut->getSkipRows() + 1))
			< trackY) {
		event.fTrackCentralRegion = event::kBad;
	}

	const DUT* dut = config.getDut(event.dut->getDUTid());
	//check if track is inside the dead region
	int col, row, pixelType;
	col = (int) ((event.trackX + 0.5 * event.dut->getPitchX())
			/ event.dut->getPitchX());
	row = (int) ((event.trackY + 0.5 * event.dut->getPitchY())
			/ event.dut->getPitchY());

	pixelType = 1;
	if (0 <= col && col <= 81 && 0 <= row && row <= 337) {
		pixelType = dut->maskedpixels[col + 1][row + 1];
	}
//  cout << "maskedpixels [" << col << "][" << row << "] =  " << pixelType << " track is = " << event.fTrackRegion << endl; 

	if (0 < pixelType && (event.fTrackRegion = event::kGood)) {

		event.fTrackRegion = event::kBad;
	}
//  cout << "TrackX =  " << trackX <<endl;
//  cout << "TrackY =  " << trackY <<endl;
//  cout << "pixelXmin =  " << pitchX*(event.dut->getSkipCols()-1) <<endl;
//  cout << "pixelXmax =  " << pitchX*(event.dut->getNcols() - (event.dut->getSkipCols()+1)) <<endl;
//  cout << "pixelYmin =  " << pitchY*(event.dut->getSkipRows()-1) <<endl;
//  cout << "pixelYmax =  " << pitchY*(event.dut->getNrows() - (event.dut->getSkipRows()+1)) <<endl;
//  cout << "event.fTrackCentralRegion =  " << event.fTrackCentralRegion <<endl;

}

