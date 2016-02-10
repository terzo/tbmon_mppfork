#ifndef SIMDATAKEEPER_H
#define SIMDATAKEEPER_H

#include <vector>

#include "simTruthHit.h"
#include "simPixelEdep.h"
//#include "tbconfig.h"

/*! \brief SIMULATION: This class keeps all "extra" event properties that exists in the case of simulated data (for one single DUT).
 *
 * Kyrre N. Sjobak (k.n.sjobak@fys.uio.no)
 */
class simDataKeeper {
private:
public:

  vector<simTruthHit*>* allTruthHits;//All truth hits in this event
  vector<simTruthHit*> truthHits;    //Truth hits in current DUT
  vector<simPixelEdep> edeps;

  simDataKeeper() {
    allTruthHits = NULL; // This will be a pointer to simTruthBuilder::truthHits
  }

  void clear() {
    truthHits.clear(); //Data is deleted by eventBuilder
    edeps.clear();
  }

 /* vector<simTruthHit*> getHitsByIden(int iden, TbConfig& config) {

    vector<simTruthHit*> retval;

    string simModName = config.simIdenNameMap_all[iden];
    for (vector<simTruthHit*>::iterator truthHit = truthHits.begin();
         truthHit != truthHits.end(); truthHit++) {
      if ( (*truthHit)->planeID == simModName ){
        retval.push_back((*truthHit));
      }
    }
    return retval;
  }
*/

};

#endif //SIMDATAKEEPER_H
