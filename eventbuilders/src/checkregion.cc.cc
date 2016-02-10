#include "checkregion.h"

void CheckRegion::buildEvent(Event &event, map<int,Event> &, TbConfig & config) {
	//Check angles in plane with iden 0
	event.fTrackMaskedRegion = event::kGood;

	bool central = true;
	double trackX(event.trackX), trackY(event.trackY);
	if( (trackX - 0.5*event.dut->getPitchX()) < (event.dut->getPitchX() * (event.dut->getSkipCols()-1)) ) central = false;
	if( (trackX - 0.5*event.dut->getPitchX()) > (event.dut->getPitchX() * (event.dut->getNcols() - (event.dut->getSkipCols()+1))) ) central = false;
	if( (trackY - 0.5*event.dut->getPitchY()) < (event.dut->getPitchY() * (event.dut->getSkipRows()-1)) ) central = false;
	if( (trackY - 0.5*event.dut->getPitchY()) > (event.dut->getPitchY() * (event.dut->getNrows() - (event.dut->getSkipRows()+1))) ) central = false;
	event.fTrackRegion = event::kGood;
	event.fTrackCentralRegion = event::kGood;
	if(not central){
		event.fTrackCentralRegion = event::kBad;
		event.fTrackRegion = event::kBad;
	}

	// (*mm).first and (*mm).second had to be switched after correcting this in
	// DUT::addMasks()
	for(vector< pair<int,int> >::iterator mm = event.dut->masks.begin();
			mm != event.dut->masks.end(); mm++){
		if( fabs(event.trackX - (event.dut->getPitchX() * (*mm).first ) ) < 1.5 * event.dut->getPitchX() and
	fabs(event.trackY - (event.dut->getPitchY() * (*mm).second )) < 1.5 * event.dut->getPitchY()) {
	 		if(config.logLevel >= kDEBUG3) { cout << "Masked track!" << endl;}
				event.fTrackMaskedRegion = event::kBad;
				event.fTrackRegion = event::kBad;
		}
	}
	if(BadRegion(event,config) == true ) {
		if(config.logLevel >= kDEBUG3) {
			cout << "Track in Bad Region!" << endl;
		}
		event.fTrackRegion = event::kBad;
	}
}

bool CheckRegion::BadRegion(Event& event, TbConfig& config) {
	double trackXT=event.trackX+0.5*event.dut->getPitchX();
	double trackYT=event.trackY+0.5*event.dut->getPitchY();
	if(strcmp(config.tbslot, "may2009") == 0 ||
			strcmp(config.tbslot, "may2009old1") == 0 ) {
		if(event.dut->getDUTid()==172) {
	 		if(trackYT>98*event.dut->getPitchY() && trackYT<142*event.dut->getPitchY())
	 			return true;    
		}
	}
	return false;
}
