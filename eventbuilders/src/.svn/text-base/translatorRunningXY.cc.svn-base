#include "translatorRunningXY.h"

void TranslatorRunningXY::initRun(TbConfig& config){
    TBALOG (kDEBUG3) << "Reading in alignment database" << std::endl;

	_dbFile = new TFile(m_filenames[config.currentRun].c_str(),"READ");
	PostShiftAlign align(_dbFile);
	for (int i = 0; i < align.getN(); i++) {
		align.getEntry(i);
		SimpleAlign *oneAlign = (SimpleAlign*) _alignColl.get(align._sensorId);
		if (!oneAlign) {
			oneAlign = new SimpleAlign(align._sensorId);
			_alignColl.add(oneAlign);
		}
		oneAlign->add(align._eventNum, align._shiftX, align._shiftY);
	}
	TBALOG (kINFO) << "Finished reading in alignment database " << m_filenames[config.currentRun].c_str() << " for run " << config.currentRun << std::endl;
}

void TranslatorRunningXY::buildEvent(Event &event, map<int,Event> &, TbConfig &config){
	float shiftX, shiftY;
	int sensorId = event.dut->getDUTid();
	SimpleAlign * myAlign = (SimpleAlign*) (_alignColl.get(sensorId));
	if (!myAlign) return;
	myAlign->getShiftXY(config.currentEntry, shiftX, shiftY);
	event.trackX = event.trackX + shiftX;
	event.trackY = event.trackY + shiftY;
}

void TranslatorRunningXY::finalizeRun(TbConfig &config){
	_alignColl.clear();
}
 
void TranslatorRunningXY::addTranslation(std::string filename, int iden, long int runNum){
	m_filenames[runNum]=filename;
}
	







