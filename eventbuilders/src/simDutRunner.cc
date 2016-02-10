#include "simDutRunner.h"

using namespace std;

void simDutRunner::init(TbConfig& config, const Event &event) {
  if (config.logLevel >= kDEBUG) {
    cout << "[ simDutRunner ]; In init()" << endl;
  }
  modelName = config.cmdLineExtras_argGetter("B_simDutRunner_model", string("passthrough"),
      "Name of model wanted. Leave to default (\"passthrough\") to use PllHits from .bdt");

  if (modelName == "passthrough") {
    dutSim = NULL;
  }
  else if (modelName == "pixel_simple") {
    dutSim = new pixel_simple(event);
  }
  else if (modelName == "Full3D_HP") {
    dutSim = new Full3D_HP(event);
  }
  else if (modelName == "Full3D_Vadim") {
    dutSim = new Full3D_Vadim(event);
  }
  else if (modelName == "passthrough") {
    dutSim = NULL;
  }
  else {
    cout << "[ simDutRunner ]; Error: Unrecognized model name \""
        << modelName << "\"" << endl;
    exit(-1);
  }


  if (dutSim != NULL) {
    dutSim->init(config);
  }

  if (config.logLevel >= kDEBUG) {
    cout << "[ simDutRunner ]; Initialized model named \""
         << modelName << "\"" << endl;
  }


}

void simDutRunner::buildEvent(Event &event, map<int,Event> &events, TbConfig &config) {
  if (dutSim == NULL) return; //Passthrough mode?
  if (config.logLevel >= kDEBUG3) {
    cout << "[ simDutRunner ]; In BuildEvent()" << endl;
  }
  event.rawHits.clear();
  dutSim->digitize(event,config,event.rawHits);

  if (config.logLevel >= kDEBUG3) {
    cout << "[ simDutRunner ]; Exiting BuildEvent()" << endl;
    for (vector<PllHit*>::iterator it = event.rawHits.begin();
        it != event.rawHits.end(); it++) {
      cout << "\t" << (*it)->col << endl;
    }
  }
}

void simDutRunner::finalize(TbConfig& config) {
  if (dutSim == NULL) return; //Passthrough mode?
  if (config.logLevel >= kDEBUG3) {
      cout << "[ simDutRunner ]; In finalize()" << endl;
  }

  dutSim->finalize(config);
}
