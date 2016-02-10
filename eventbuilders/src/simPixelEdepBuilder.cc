#include "simPixelEdepBuilder.h"

using namespace std;

void simPixelEdepBuilder::initRun(TbConfig &config) {
  if (config.logLevel >= kDEBUG3) {
      cout << "[ simPixelEdepBuilder ]; In initRun()" << endl;
  }
  //Safety
  if (not config.isSimulation) {
    cout << "[ simPixelEdepBuilder ]; Error: This is not simulation data!" << endl;
    exit(-1);
  }
  if (strcmp(config.simDataPath_subdir, "notset_inBuilder") == 0) {
    cout << "[ simPixelEdepBuilder ]; Error: simBaseBuilder has to run first!" << endl;
    exit(-1);
  }
  /*
  if (not config.fNeedSimPixelEdep) {
    if (config.logLevel >= kDEBUG) {
      cout << "[ simPixelEdepBuilder ]; config.fNeedSimPixelEdep == false, skipping builder" << endl;
    }
    return;
  }
  */

  //Open files, get folder name from config
  char* fileName = new char[600];
  for (map<int,string>::iterator it = config.simIdenNameMap.begin();
      it != config.simIdenNameMap.end(); it++) {

    sprintf(fileName, "%sTBsim_moduleDUT_%s.dat", config.simDataPath_subdir, (*it).second.c_str());
    if (config.logLevel >= kDEBUG) {
        cout << "[ simPixelEdepBuilder ]; Opening edepfile " << fileName << endl;
    }
    pixelFiles[(*it).first] = new ifstream(fileName);


    if (! pixelFiles[(*it).first]->is_open()) {
      cout << "[ simPixelEdepBuilder ]; Error: edepfile named " << fileName << " could NOT be opened!" << endl;
      exit(-1);
    }

    //Create edep vector for this iden;
    //edeps[(*it).first] = vector<simPixelEdep>();
  }
  delete fileName;

  if (config.logLevel >= kDEBUG3) {
      cout << "[ simPixelEdepBuilder ]; Exiting initRun()" << endl;
  }
}

void simPixelEdepBuilder::initEvent(TbConfig &config, map<int,Event> &events){

};

void simPixelEdepBuilder::buildEvent(Event &event, map<int,Event> &events, TbConfig & config){
  if (config.logLevel >= kDEBUG3) {
      cout << "[ simPixelEdepBuilder ]; In buildEvent()" << endl;
  }

  //if (not config.fNeedSimEdep) return;

  //Prepare storage
  int iden = event.dut->getDUTid();
  ifstream* edepFile = pixelFiles[iden];
  vector<simPixelEdep>& edepVector = event.simData->edeps;

  //Fetch data
  int eventNum(-1), numEntries(-1);
  bool seekEvent(true);
  do {
    //Fuglyness: Need to use dereference operator on ifstream
    *edepFile >> eventNum; *edepFile >> numEntries; //Entry header
    if (eventNum == config.simCurrentEventNum) {
      //Found correct simulation event, parse it
      seekEvent = false;
      for (int i = 0; i < numEntries && edepFile->good(); i++) {
        simPixelEdep edep = simPixelEdep();
        *edepFile >> edep.edep;
        *edepFile >> edep.pos;
        *edepFile >> edep.posLocal;
        edep.posLocal.localTranslate(event);
        //cout << eventNum << "\t" << config.currentEntry << endl;
        //edep.print();
        edepVector.push_back(edep);
      }
      if (not edepFile->good()) {
        cout <<"[ simPixelEdepBuilder ]; Error: Bad file while parsing event" << endl;
        exit(-1);
      }
    }
    else {
      //Skip this event. +1 b.c. first line got is just line ending
      // from header line
      for (int i = 0; i < numEntries+1 && edepFile->good(); i++) {
        edepFile->getline(bitBucket,256);
      }
      if (not edepFile->good()) {
        cout <<"[ simPixelEdepBuilder ]; Error: Bad file while filling bitBucket" << endl;
        exit(-1);
      }
    }
  } while (seekEvent && edepFile->good());
  if (not edepFile->good()) {
    cout << "[ simPixelEdepBuilder ]; Error: Bad file while searching edepFile" << endl;
    exit(-1);
  }

  if (config.logLevel >= kDEBUG3) {
      cout << "[ simPixelEdepBuilder ]; Exiting buildEvent()" << endl;
  }
}

void simPixelEdepBuilder::finalizeRun(TbConfig& config) {
  //if (not config.fNeedSimEdep) return;

  //Close files, delete pointers
  if (config.logLevel >= kDEBUG) {
    cout << "[ simEdepBuilder ]; Closing edepfiles" << endl;
  }

  for (map<int,ifstream*>::iterator it = pixelFiles.begin();
      it != pixelFiles.end(); it++) {
    (*it).second->close();
    delete (*it).second;
  }
  pixelFiles.clear();
}
