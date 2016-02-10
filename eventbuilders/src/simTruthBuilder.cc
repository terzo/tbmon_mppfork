#include "simTruthBuilder.h"

using namespace std;

void simTruthBuilder::init(TbConfig& config) {
  useTruthTracking = config.cmdLineExtras_argGetter("B_simTruthBuilder_truthTrack", false, "Use truth info for position estimation in devices");
}

void simTruthBuilder::initRun(TbConfig& config) {
  //Safety
  if (not config.isSimulation) {
    cout << "[ simTruthBuilder ]; Error: This is not simulation data!" << endl;
    exit(-1);
  }
  if (strcmp(config.simDataPath_subdir, "notset_inBuilder") == 0) {
    cout << "[ simTruthBuilder ]; Error: simBaseBuilder has to run first!" << endl;
    exit(-1);
  }
  /*
  if (not config.fNeedSimTruth) {
    if (config.logLevel >= kDEBUG) {
      cout << "[ simPixelEdepBuilder ]; config.fNeedSimTruth == false, skipping builder" << endl;
    }
    return;
  }
  */

  //Open file, get folder name from config
  char* fileName = new char[600];
  sprintf(fileName, "%sTBsim_truth.dat", config.simDataPath_subdir);
  if (config.logLevel >= kDEBUG) {
    cout << "[ simTruthBuilder ]; Opening truthFile " << fileName << endl;
  }
  truthFile.open(fileName, ios::in);
  if (!truthFile.is_open()) {
    cout << "[ simTruthBuilder ]; Error: truthFile named " << fileName << " could NOT be opened!" << endl;
    exit(-1);
  }
  delete fileName;
}

void simTruthBuilder::initEvent(TbConfig &config, const Event &event, map<int,Event> &events) {
  if (config.logLevel >= kDEBUG3) {
      cout << "[ simTruthBuilder ]; In initEvent()" << endl;
  }
  //if (not config.fNeedSimEdep) return;


  //Clear any data from previous run
  for (vector<simTruthHit*>::iterator it = truthHits.begin(); it != truthHits.end(); it++) {
    delete (*it);
  }
  truthHits.clear();

  // Fetch the truth data for the correct event from file,
  // store it in this->truthHits
  int eventNum(-1), numEntries(-1);
  bool seekEvent(true);
  do {
    truthFile >> eventNum; truthFile >> numEntries; //Entry header
    if (eventNum == config.simCurrentEventNum) {
      //Found correct simulation event, parse it
      seekEvent = false;
      for (int i = 0; i < numEntries && truthFile.good(); i++) {
        simTruthHit* hit = new simTruthHit();
        truthFile >> hit->pos;
        truthFile >> hit->posLocalRaw;
        hit->posLocal = hit->posLocalRaw;
        hit->posLocal.localTranslate(event);
        truthFile >> hit->momDir;
        truthFile >> hit->momDirLocal;
        truthFile >> hit->kinE;
        truthFile >> hit->planeID;
        truthFile >> hit->firstStep;
        truthFile >> hit->particleID;
        truthFile >> hit->particleType;
        truthFile >> hit->stepNum;
        //hit->print();
        truthHits.push_back(hit);
      }
      if (not truthFile.good()) {
        cout <<"[ simTruthBuilder ]; Error: Bad file while parsing event" << endl;
        exit(-1);
      }
    }
    else {
      //Skip this event. +1 b.c. first line got is just line ending
      // from header line
      for (int i = 0; i < numEntries+1 && truthFile.good(); i++) {
        truthFile.getline(bitBucket,256);
      }
      if (not truthFile.good()) {
        cout <<"[ simTruthBuilder ]; Error: Bad file while filling bitBucket" << endl;
        exit(-1);
      }
    }
  } while (seekEvent && truthFile.good());
  if (not truthFile.good()) {
    cout << "[ simTruthBuilder ]; Error: Bad file while searching truthFile" << endl;
    exit(-1);
  }

  if (config.logLevel >= kDEBUG3) {
        cout << "[ simTruthBuilder ]; Exiting initEvent()" << endl;
  }
}



void simTruthBuilder::buildEvent(Event &event, map<int,Event> &, TbConfig & config) {
  //if (not config.fNeedSimEdep) return;

  //Needed for accessing telescope truth huts
  event.simData->allTruthHits = &(this->truthHits);

  // Push the truth data from this->truthHits onto the event.
  // Use some mapping between DUT iden and TestBeamSim naming loaded from config and
  // specified in builder
  string simModName = config.simIdenNameMap[event.dut->getDUTid()];
  for (vector<simTruthHit*>::iterator truthHit = truthHits.begin();
       truthHit != truthHits.end(); truthHit++) {
    if ( (*truthHit)->planeID == simModName ){
      event.simData->truthHits.push_back((*truthHit));
    }
  }

  if(useTruthTracking) {
    if(event.simData->truthHits.size() != 1) {
      event.fTrack = event::kBad;
      return;
    }
    simTruthHit* hit = event.simData->truthHits[0];
    //Sim (global) X = pixel Y and the other way around
    double dz = hit->momDirLocal.data[2];
    event.setEvent(hit->posLocal.data[0],hit->posLocal.data[1],
        hit->momDirLocal.data[1]/dz,hit->momDirLocal.data[0]/dz,0,0);
  }
}


void simTruthBuilder::finalizeRun(TbConfig& config) {
  //if (not config.fNeedSimEdep) return;

  //Close file
  if (config.logLevel >= kDEBUG) {
    cout << "[ simTruthBuilder ]; Closing truthFile" << endl;
  }
  truthFile.close();
}
