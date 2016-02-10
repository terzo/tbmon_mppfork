#include "totcalibreader.h"

void ToTCalibReader::init(TbConfig& config) {
	//char key [100], description[100];
	//cout << config.currentRun << endl;
	//sprintf(key,"B_totcalibreader_11193",config.currentRun);
	//sprintf(description,"Calibration file path run %u",config.currentRun);
	string calibpath = config.cmdLineExtras_argGetter("B_totcalibreader",
			string(""), "Calibration file path");
	if (calibpath == "")
		return;

	for (list<DUT*>::const_iterator i = config.dutList.begin();
			i != config.dutList.end(); i++) {
		(*i)->addToTCalib(calibpath);
	}

}
