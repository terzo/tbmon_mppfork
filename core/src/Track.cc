#include "Track.h"

ClassImp(BatCluster)
ClassImp(PllHit)
ClassImp(HptHit)
ClassImp(TrackParams)
ClassImp(TrackErrors)
ClassImp(Track)

Track::Track() :
  iden(0),
  trig(0),
  nBatCluster(0),  batCluster(TClonesArray("BatCluster")),
  nPllHit(0),      pllHit(TClonesArray("PllHit")),
  nHptHit(0),      hptHit(TClonesArray("HptHit")),
  nTrackParams(0), trackParams(TClonesArray("TrackParams")),
  nTrackErrors(0), trackErrors(TClonesArray("TrackErrors"))
{;}
