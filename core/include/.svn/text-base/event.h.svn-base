#ifndef EVENT_H
#define EVENT_H

#include "Track.h"
#include "dut.h"
#include "simDataKeeper.h"
//#include <list> //Not used
#include <vector>
#include <iostream>

using namespace std;

namespace event {
typedef vector<PllHit*> cluster_t;
enum eventflag_t {
	kBad, kGood, kUnknown
};

}
;

/*! \brief Contains data of one event for one DUT (respectively links):
 *
 */
class Event {
public:
	//Base event
	DUT* dut;
	Track* track;
	vector<PllHit*> hits;
	vector<PllHit*> rawHits;
	double trackX;
	double trackY;
	double dxdz;
	double dydz;
	double chi2;
	double ndof;

	//global parameters
	float gr_RotXY, gr_RotZX, gr_RotZY;

	//Clusters
	vector<event::cluster_t> clusters;
	vector<event::cluster_t> rawClusters;
	vector<event::cluster_t> euClusters;

	//Simulation
	simDataKeeper* simData;

	//Flags
	/** A list of all the flags pertaining to an Event object and the (rough) code lines where they are potentially set. A block of several close lines may be referenced only once. Bracketed () references mean not included in the current list of eventbuilders in config. */
	event::eventflag_t fBase; ///< event.h 77; (BAT 39) - guaranteed kGood if setEvent called
	event::eventflag_t fClusters; // Are clusters built? ///< CLUSTF 23 - guaranteed kGood in buildEvent of ClusterFinder
	event::eventflag_t fTrack; ///< event.h 78; EUBT 160; (CHI2 13); (BAT 36); (ANGLE 55,59) - by default kGood, unless fewer matched hits are found than m_nMatch EuBuildTrack? param (set in our driver.cc as 1)
	event::eventflag_t fTrackChi2; // Is track Chi2 accepted? ///< (CHI2 12)
	event::eventflag_t fTrackAngle; ///< (ANGLE 43,51,54,58)
	event::eventflag_t fHits; ///< completely redundant, neither set nor used
	event::eventflag_t fTrackCentralRegion; ///< CREG 18,20 - kBad if track passes through pixel inside border specified by event::skipCols, event::skipRows (hard coded as 2 and 16 respectively)
	event::eventflag_t fTrackMaskedRegion; ///< CREG 9,29 - kBad if the track matches (within 1.5 * pitch dimensions) a pixel which is masked out by DUT mask
	event::eventflag_t fTrackRegion; ///< CREG 21,30,35 - kBad if fTrackMaskedRegion or fTrackCentralRegion kBad
	event::eventflag_t fTrackMatchNeighbour; ///< EUBT 160 - same as fTrack
	event::eventflag_t fEtaCorrections; ///< CLUSTF 24
	event::eventflag_t fDutSync; ///< (DUTSYNC 56,59,64)
	event::eventflag_t fSimSync; //Is simulation and reconstructed data well synced up? ///< (SIMBASE 81,98)
	event::eventflag_t fEtaCut; //BAT eta cuts ///< (ETACUT 43 *******due to a presumed error line 46, the flag is never set kGood********)

	Event(DUT* dut = NULL) :
			dut(dut), track(NULL), simData(NULL) {
	}
	;

	void setEvent(double x, double y, double dxdz, double dydz, double chi2,
			int ndof) {
		this->trackX = x;
		this->trackY = y;
		this->dxdz = dxdz;
		this->dydz = dydz;
		this->chi2 = chi2;
		this->ndof = ndof;
		fBase = event::kGood;
		fTrack = event::kGood;
	}
	;
	void clear() {
		track = NULL;
		hits.clear();
		rawHits.clear();
		clusters.clear();
		rawClusters.clear();
		if (simData != NULL)
			simData->clear();
		fBase = event::kUnknown;
		fClusters = event::kUnknown;
		fTrackChi2 = event::kUnknown;
		fTrackAngle = event::kUnknown;
		fTrack = event::kUnknown;
		fTrackCentralRegion = event::kUnknown;
		fTrackRegion = event::kUnknown;
		fTrackMaskedRegion = event::kUnknown;
		fTrackMatchNeighbour = event::kUnknown;
		fEtaCorrections = event::kUnknown;
		fDutSync = event::kUnknown;
		fSimSync = event::kUnknown;
		fEtaCut = event::kUnknown;
		chi2 = 0;
		ndof = 0;
	}
	;
};
#endif //EVENT_H
