#ifndef LOOPER_H
#define LOOPER_H

// standard header files
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <vector>

// root header files
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

class Looper;

// tbmon header files
#include "tbconfig.h"
#include "Track.h"
#include "event.h"
#include "tbanalysis.h"
#include "dut.h"

using namespace std;

/*! \brief Manages the whole event loop by calling event builders and analysis processors. Is manged on
 * its own by tbconfig.
 *
 */
class Looper {
public:
	// main loop over all runs, calls initRun, eventLoop, finalizeRun for each run
	// and after the loop finalize
	void loop(TbConfig &config);
	// loop over event in a run, calls event
	void eventLoop(TbConfig &config);
	// virtual functions
	void event(TbConfig &config);
	void finalize(TbConfig &config);
	void finalizeRun(TbConfig &config);
	// calls initRun for alls DUTs, EventBuilders and Analyses
	void initRun(TbConfig &config);
};

#endif //LOOPER_H
