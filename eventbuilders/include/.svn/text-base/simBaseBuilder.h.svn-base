#ifndef SIMBASEBUILDER_H
#define SIMBASEBUILDER_H

#include <dirent.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "battrack.h"
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

using namespace std;

/*! \brief SIMULATION: This event builder sets up basic configurations for all simulation-specific eventBuilders. It has to run before other simulation-specific builders!
 *
 * Does work that for "real data" is done directly in Looper::loop()
 *
 * Kyrre N. Sjobak
 */
class simBaseBuilder : public EventBuilder {
public:
  simBaseBuilder() {
    name = "simBaseBuilder";
  };
  virtual void initRun(TbConfig &config);
  virtual void initEvent(TbConfig &config, map<int,Event> &events); //Stuff that is common to all DUTs
  virtual void buildEvent(Event &event, map<int,Event> &events, TbConfig &); //Per-DUT building
  virtual void finalizeRun(TbConfig& config);

private:
  ifstream syncMapFile;
};



#endif //SIMBASEBUILDER_H
