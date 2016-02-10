#include "simBaseBuilder.h"

using namespace std;


void simBaseBuilder::initRun(TbConfig& config) {
  //Safety
  if (not config.isSimulation) {
    cout << "[ simBaseBuilder ]; Error: This is not simulation data!" << endl;
    exit(-1);
  }

  //Figure out name of subfolder where data is stored, save to config
  DIR* d = NULL;
  d = opendir(config.simDataPath);
  if (d == NULL) {
    cout << "[ simBaseBuilder ]; Error: Could not open directory" << config.simDataPath << endl;
    exit(-1);
  }
  struct dirent* de = NULL;
  char* currentRun_baseName = new char [100];
  sprintf(currentRun_baseName, "run%u_", config.currentRun);
  while(de = readdir(d)) {
    if (de->d_type != DT_DIR or
        strncmp(de->d_name, currentRun_baseName, strlen(currentRun_baseName)) != 0)
      continue;

    sprintf(config.simDataPath_subdir, "%s%s/", config.simDataPath, de->d_name);

    if (config.logLevel >= kDEBUG) {
      cout << "[ simBaseBuilder ]; Found simulation subdir " << de->d_name << endl;
      cout << "[ simBaseBuilder ]; => simulation datadir  " << config.simDataPath_subdir << endl;
    }

  }
  closedir(d);

  //Get TestBeamSim<->BDT event number map digiSyncMap.dat
  char* fileName = new char[600];
  sprintf(fileName, "%sdigiSyncMap.dat", config.simDataPath_subdir);
  if (config.logLevel >= kDEBUG) {
    cout << "[ simBaseBuilder ]; Opening digiSyncMap " << fileName << endl;
  }
  syncMapFile.open(fileName, ios::in);
  if (!syncMapFile.is_open()) {
    cout << "[ simTruthBuilder ]; Error: digiSyncMap file named " << fileName << " could NOT be opened!" << endl;
    exit(-1);
  }
  delete fileName;

}
void simBaseBuilder::initEvent(TbConfig &config, map<int,Event> &events) {
  if (config.logLevel >= kDEBUG3) {
      cout << "[ simBaseBuilder ]; In initEvent()" << endl;
  }
  //Get BatTrack's Track object
  Track* m_track = ((BatTrack*) config.eventbuilders_map["BatTrack"])->getTrack();

  //Check sync
  int evNum = m_track->trig; //Get main event number
  //Check that BAT's are in sync with main event number
  for (int i = 0; i < m_track->nBatCluster; i++) {
    int thisEvNum = ((BatCluster*) m_track->batCluster[i])->trig;
    if (evNum != thisEvNum) {
      cout << "[ simBaseBuilder ]; Error: Found disagreeing trigger numbers in BAT clusters" << endl;
      exit(-1);
    }
  }
  if (evNum == -1) {
    //No BAT cluster found, skip this event
    for (map<int,Event>::iterator it = events.begin(); it != events.end(); it++) {
      (*it).second.fSimSync = event::kBad;
    }
    if (config.logLevel >= kDEBUG3) {
      cout << "[ simBaseBuilder ]; Bad simSync (no batClusters)" << endl;
    }
    return;
  }
  //Check that DUTs are in sync with main event number
  for (int i = 0; i < m_track->nPllHit; i++) {
    int thisEvNum = ((PllHit*) m_track->pllHit[i])->trig;
    if (evNum != thisEvNum) {
      cout << "[ simBaseBuilder ]; Error: Found disagreeing trigger number in PllHit" << endl;
      exit(-1);
    }
  }
  //If got to this point, everything is good :)
  for (map<int,Event>::iterator it = events.begin(); it != events.end(); it++) {
        (*it).second.fSimSync = event::kGood;
  }

  //Read synchronization info from file
  int digiTrig, simTrig;
  do {
    syncMapFile >> digiTrig;
    syncMapFile >> simTrig;
  } while(digiTrig != evNum);
  config.simCurrentEventNum = simTrig;

  if (config.logLevel >= kDEBUG3) {
      cout << "[ simBaseBuilder ]; Exiting initEvent(); digiTrig="
           << digiTrig << " simTrig=" << simTrig << endl;
  }

}
void simBaseBuilder::buildEvent(Event &event, map<int,Event> &events, TbConfig &) {
  if (event.simData == NULL) {
    event.simData = new simDataKeeper();
  }
}

void simBaseBuilder::finalizeRun(TbConfig& config) {
  //Close file
  if (config.logLevel >= kDEBUG) {
    cout << "[ simBaseBuilder ]; Closing digiSyncMap file" << endl;
  }
  syncMapFile.close();
}
