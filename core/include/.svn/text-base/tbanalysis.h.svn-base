#ifndef TBANALYSIS_H
#define TBANALYSIS_H
#define TBALOG(level) if(level <= config.logLevel) cout << "[ " << this->name << " ]: " 

class TbAnalysis;
#include "event.h"
#include "tbconfig.h"
#include "dut.h"
#include <map>
#include <iostream>
#include "tbutils.h"
#include <sstream>

using namespace std;

/*! \brief Mostly virtual base class of all analysis processors.
 *
 */
class TbAnalysis {
public:
	char* basename;
	char* name;
	int iden;
	void setDut(const DUT* dut) {
		this->iden = dut->getDUTid();
	}
	void setName(const char* name) {
		this->basename = (char*) name;
		this->name = new char[500];
		sprintf(this->name, "%s-%i", name, iden);
	}
	TbAnalysis() {
		name = NULL;
		basename = NULL;
	}
	;
	//Stuff to be done before anything else
	//virtual void init(TbConfig &config) = 0;
	virtual void init(TbConfig &config) = 0;
	//How to process an event
	virtual void event(const TbConfig &config, const Event &event) = 0;
	//Stuff to be done after all events have been processed
	//virtual void finalize(const TbConfig &config) = 0;
	virtual void finalize(const TbConfig &config) = 0;
	//Stuff to be done at the beginning of each run
	virtual void initRun(const TbConfig &config) {
		;
	}
	//Stuff to be done at the end of each run
	virtual void finalizeRun(const TbConfig &config) {
		;
	}
};

#endif //TBANALYSIS_H
