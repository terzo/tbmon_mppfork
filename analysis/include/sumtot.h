#ifndef SUMTOT_H
#define SUMTOT_H

#include "event.h"
#include "tbanalysis.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TVector2.h"
#include <string>

class SumTot : public TbAnalysis{

private:

  TH1D *h_fullChipTot;
  TH1D *h_fullChipQ;
  TH1D* h_maxClusterTot;
  TH1D* h_maxClusterQ;
  TH1D *h_maxCellTot;
  TH1D *h_maxCellQ;
  TH1D* h_matchClusterTot;
  TH1D* h_matchPixelTot;
  TH1D* h_matchClusterQ;
  TH1D* h_matchClusterQ_1;
  TH1D* h_matchClusterQ_2;
  TH1D* h_matchClusterQ_3;
  TH1D* h_matchClusterQ_4;
  TH1D* h_matchClusterQ_5plus;
  TH1D* h_clusterq_no_bg;
  TH2D *h_elecMap;  
  TH3D *h_elecMapTot;
  TH1D* h_elecDist;
  TH2D* h_elecDistTot;

  TH3D* h_matchClusterTotMap;

  TH3D* h_matchClusterQMap;
  TH2D *h_mapElec2;  
  
  std::map< int, std::map<int,TH1D*> > hh_matchPixelTotCol, hh_pureMatchPixelTotCol, hh_shareMatchPixelTotCol, hh_matchPixelTotRow, hh_pureMatchPixelTotRow;



  bool doCharge;

  void totMap( const TbConfig& config, const char* analysisName, const char* histoName,TH3D* totmap3d ) const;
  void ElecDistTot( const TbConfig& config, const char* analysisName, const char* histoName, TH2D* totelecdist ) const;
  void totElecMap( const TbConfig& config, const char* analysisName, const char* histoName, TH3D* totmap3d ) const;

public:  
  virtual void init(TbConfig &config);
  virtual void event(const TbConfig &config, const Event &event);
  virtual void finalize(const TbConfig &config); 

};
#endif
