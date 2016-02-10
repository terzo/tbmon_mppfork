#include "translator.h"

void Translator::initRun(TbConfig& config) {
	m_translations.clear();
	for (list<TransCalib>::iterator it = m_calibs.begin(); it != m_calibs.end();
			it++) {
		if (config.currentRun >= (*it).firstRun
				and config.currentRun <= (*it).lastRun) {
			m_translations[(*it).iden] = (*it);
		}
	}
}

void Translator::buildEvent(Event &event, map<int, Event> &, TbConfig &) {
	map<int, TransCalib>::iterator translations = m_translations.find(
			event.dut->getDUTid());
	if (translations == m_translations.end()) {
		return;
	}
	event.trackX = event.trackX - (*translations).second.shiftX;
	event.trackY = event.trackY - (*translations).second.shiftY;
}

void Translator::addTranslation(int iden, int firstRun, int lastRun,
		double shiftX, double shiftY) {
	TransCalib calib;
	calib.iden = iden;
	calib.firstRun = firstRun;
	calib.lastRun = lastRun;
	calib.shiftX = shiftX;
	calib.shiftY = shiftY;
	m_calibs.push_back(calib);
}

void Translator::addTranslation(const char* filename, int iden, double addX,
		double addY) {
	ifstream file(filename);
	int runnumber;
	while (!file.eof()) {
		TransCalib calib;
		file >> runnumber >> calib.shiftX >> calib.shiftY;
		calib.shiftX += addX;
		calib.shiftY += addY;
		calib.iden = iden;
		calib.firstRun = runnumber;
		calib.lastRun = runnumber;
		m_calibs.push_back(calib);
	}
	file.close();
}

