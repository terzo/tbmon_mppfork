#ifndef TBTRACK_H
#define TBTRACK_H

#include "TClonesArray.h"

/*! \brief Contains one cluster for BAT data
 *
 */
class BatCluster : public TObject {
public:
  BatCluster() {;}
  Int_t trig; // Trigger number
  Int_t iden; // Module identifier
  Int_t claye; // Layer: 0 vertical 1 horizontal
  Float_t cstrp; // Cluster barycenter
  Float_t ccnts; // Summed ADC counts
  Float_t cnois; // Summed strip noise
  Float_t ccorr; // Eta-corrected barycenter

  ClassDef(BatCluster, 1)
};
/*! \brief Contains one DUT hit
 *
 */
class PllHit : public TObject {
public:
  PllHit() {;}
  Int_t trig; // Trigger number
  Int_t iden; // Module identifier
  Int_t chp;  // Chip number
  Int_t row;  // Pixel row (0-159)
  Int_t col;  // Pixel column (0-17)
  Int_t tot;  // Channel ToT (0-255)
  Int_t lv1;  // BX id relative to lv1 (0-15)

  ClassDef(PllHit, 1)
};
/*! \brief Contains probably information of one hit for BAT telescope planes(?) It is not used anywhere,
 * candiate for deletion.
 *
 */
class HptHit : public TObject {
public:
  HptHit() {;}
  Int_t trig; // Trigger number
  Int_t iden; // Module identifier
  Int_t chn;  // Channel number
  Int_t hit;  // TDC count

  ClassDef(HptHit, 1)
};
/*! \brief Contains track information for BAT data.
 *
 */
class TrackParams : public TObject {
public:
  TrackParams() {;}
  Int_t iden;         // Reference plane identifier
  #ifdef FOUR_PARAMS
  Float_t params[4];  // Track parameters: x0 y0 dx dy
  #endif
  #ifndef FOUR_PARAMS
  Float_t params[5];  // Track parameters: x0 y0 dx dy -|B|q/p
  #endif
  ClassDef(TrackParams, 1)
};
/*! \brief Track errors - not used anywhere, probably BAT as well?
 *
 */
class TrackErrors : public TObject {
public:
  TrackErrors() {;}
  Int_t iden;         // Reference plane identifier
  #ifdef FOUR_PARAMS
  Float_t errors[10]; // Covariance matrix: Upper triangular columns
  #endif
  #ifndef FOUR_PARAMS
  Float_t errors[15]; // Covariance matrix: Upper triangular columns
  #endif


  ClassDef(TrackErrors, 1)
};
/*! \brief Contains data of one track for one DUT
 *
 */
class Track : public TObject {
public:
  Track();
  Int_t trig; // Trigger number
  Int_t iden; // Module identifier
  Float_t ct; // Chi2 (raw)
  Float_t ft; // Hit count (weighted)
  Float_t gt; // Chi2 (weighted)
  Int_t nBatCluster;  TClonesArray batCluster;
  Int_t nPllHit;      TClonesArray pllHit;
  Int_t nHptHit;      TClonesArray hptHit;
  Int_t nTrackParams; TClonesArray trackParams;
  Int_t nTrackErrors; TClonesArray trackErrors;
  
  ClassDef(Track, 1)
};
#endif //TBTRACK_H
