#include "maskreader.h"

void MaskReader::init(TbConfig& config) {
	char key[100], description[100];
	//cout << config.currentRun << endl;
	string calibpath;

	for (list<DUT*>::const_iterator i = config.dutList.begin();
			i != config.dutList.end(); i++) {
		sprintf(key, "B_maskreader_%u", (*i)->getDUTid());
		sprintf(description, "Maskfile path DUT iden %u", (*i)->getDUTid());
		calibpath = config.cmdLineExtras_argGetter(key, string(""),
				description);
		if (calibpath == "")
			continue;
		(*i)->addMasks(calibpath.c_str());
		//cout << "Job done" << endl;		
	}

}
