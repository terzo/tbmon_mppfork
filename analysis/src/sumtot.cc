#include "sumtot.h"
#include <iostream>
#include "clusters.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"

void SumTot::init(TbConfig &config){

  const DUT* dut = config.getDut(this->iden);
  double pitchX = dut->getPitchX();
  double pitchY = dut->getPitchY();
  int nCols = dut->getNcols();
  int nRows = dut->getNrows();
  // hack, but gives reasonable range for both FE-I3 and FE-I4
  int maxToT = 250;//(int) floor(2304/nCols);

  h_fullChipTot = new TH1D("", ";Pixel hit charge [ToT]", maxToT, 0, maxToT);
  h_fullChipQ = new TH1D("", ";Pixel hit charge [e^{-}]", 150, 0, 100000);
  h_maxClusterTot = new TH1D("", ";Max cluster charge [ToT]", maxToT, 0, maxToT);
  h_maxClusterQ = new TH1D("", ";Max cluster charge [e^{-}]", 150, 0, 100000);
  h_maxCellTot = new TH1D("", ";Max cell charge [ToT]", maxToT, 0, maxToT);
  h_maxCellQ = new TH1D("", ";Max cell charge [e^{-}]", 150, 0, 100000);

  h_matchClusterTot = new TH1D("", ";Matching cluster charge [ToT]", maxToT, 0, maxToT);
  h_matchPixelTot = new TH1D("", ";Matching pixel charge [ToT]", maxToT, 0, maxToT);
  h_matchClusterQ = new TH1D("", ";Matching cluster charge [e^{-}]", 300, 0, 100000);
  h_matchClusterQ_1 = new TH1D("", ";Cluster charge [e^{-}]", 300, 0, 100000);
  h_matchClusterQ_2 = new TH1D("", ";Cluster charge [e^{-}]", 300, 0, 100000);
  h_matchClusterQ_3 = new TH1D("", ";Cluster charge [e^{-}]", 300, 0, 100000);
  h_matchClusterQ_4 = new TH1D("", ";Cluster charge [e^{-}]", 300, 0, 100000);
  h_matchClusterQ_5plus = new TH1D("", ";Cluster charge [e^{-}]", 300, 0, 100000);
  h_clusterq_no_bg = new TH1D("", ";Cluster charge [e^{-}] 0-300#mum", 300, 0, 100000);

  h_matchClusterTotMap = new TH3D("", ";Column;Row;Cluster charge [ToT]", nCols, 0, nCols, nRows, 0, nRows, maxToT, 0, maxToT);
  h_matchClusterQMap = new TH3D("", ";Column;Row;Cluster charge [e^{-}]", nCols, 0, nCols, nRows, 0, nRows, 150, 0, 100000);

  h_elecMap = new TH2D("", ";X-Distance to electrode [#mum];Y-Distance to electrode [#mum]", 50, -100, 100, pitchY, -pitchY, pitchY);
  h_elecMapTot = new TH3D("", ";X-Distance to electrode [#mum];Y-Distance to electrode [#mum];Matching cluster charge [ToT]", 50, -100., 100., pitchY, -pitchY, pitchY, 0.5*maxToT, 0, maxToT);
  h_elecDist = new TH1D("", ";Distance to center of electrode [#mum]", 50, 0, 100);
  h_elecDistTot = new TH2D("", ";Distance to center of electrode [#mum];Matching cluster charge [ToT]", 50, 0, 100., 0.5*maxToT, 0, maxToT);

  doCharge = false;
  gStyle->SetOptStat(111111);

}


void SumTot::event(const TbConfig &config, const Event &event){

  // Check that track has passed all cuts
  if( event.fTrack != event::kGood ) { return; }
  // Check that clusters have been made successfully
  if( event.fClusters != event::kGood ) { return; }
  // Require that we have hits
  if( event.hits.size() == 0 ) { return; }

  h_fullChipTot->Fill( cluster::getSumTot( event.hits ) );
  if( doCharge==true ) { h_fullChipQ->Fill( cluster::getSumCharge( event.hits, event ) ); }
  
  int maxCluster = cluster::getMaxTotCluster(event.clusters);
  if( maxCluster > -1 and cluster::getSumTot(event.clusters.at(maxCluster)) > 0 ) {
    if( cluster::isCentralRegion(event.clusters.at(maxCluster), event) ) {  
      h_maxClusterTot->Fill( cluster::getSumTot( event.clusters.at(maxCluster) ));
      if( doCharge==true ) { h_maxClusterQ->Fill( cluster::getSumCharge( event.clusters.at(maxCluster), event )); }
    }
  }

  int maxCell = cluster::getMaxTotCell( event.hits );
  if( maxCell > -1 ){
    h_maxCellTot->Fill( event.hits.at(maxCell)->tot );
    if( doCharge==true ) { h_maxCellQ->Fill( event.dut->q( event.hits.at(maxCell)->tot, event.hits.at(maxCell)->col, event.hits.at(maxCell)->row ) ); }
  }
  
  // Use only clusters on tracks in the central region and exclude bad regions
  // Are we in the central region on the chip?
  if( event.fTrackCentralRegion != event::kGood) { return; }
  TBALOG(kDEBUG2) << "Central track found!" << endl; 
  // Are we in a good region of the chip?
  if( event.fTrackRegion != event::kGood) { 
    TBALOG(kDEBUG2) << "Track in bad region, col = " << tbutils::getCol(event.trackX,event) << " row = "<< tbutils::getRow(event.trackY,event) << endl; 
    return;
  }
  TBALOG(kDEBUG2) << "Track in good region found !" << endl; 
  
  int matchCluster = cluster::getMatched( event.clusters, event );
  if( matchCluster == -1 ) {
    TBALOG(kDEBUG2) << "No matching cluster found!" << endl; 
    return;
  }
  // Make sure that all hits are in the central region
  bool centralHits = cluster::isCentralRegion( event.clusters.at(matchCluster), event );
  if( centralHits == false ) {
    TBALOG(kDEBUG2) << "Found cluster with hits outside central region!" << endl; 
    return;
  }
  // Make sure that matched cluster actually has ToT
  if( cluster::getSumTot( event.clusters.at(matchCluster) ) <= 0 ) {
    TBALOG(kDEBUG2) << "Matched cluster has ToT = "<< cluster::getSumTot( event.clusters.at(matchCluster) ) << ", skip event!" << endl; 
    return;
  }  
  TBALOG( kDEBUG2 ) << "Found matching cluster!" << endl;

  h_matchClusterTot->Fill( cluster::getSumTot( event.clusters.at(matchCluster) ));
 
  std::vector<std::pair<int,int> > posColToT, posRowToT, posSharingCol;
  bool pureRow=true, pureCol=true;
  int roww=0,colw=0,clusterCol=-1000;
  for(int ii = 0; ii < event.clusters.at(matchCluster).size(); ii++)
  {
     h_matchPixelTot->Fill( event.clusters.at(matchCluster).at(ii)->tot );
     if(ii==0) 
     {
        posColToT.push_back( std::make_pair(event.clusters.at(matchCluster).at(ii)->col, event.clusters.at(matchCluster).at(ii)->tot) );
	roww=event.clusters.at(matchCluster).at(ii)->row;
	posRowToT.push_back( std::make_pair(event.clusters.at(matchCluster).at(ii)->row, event.clusters.at(matchCluster).at(ii)->tot) );
	colw=event.clusters.at(matchCluster).at(ii)->col;
     }
     else 
     {
     	bool found=false;
        for(std::vector<std::pair<int,int> >::iterator it=posColToT.begin(); it!=posColToT.end(); ++it)
	{
           if( event.clusters.at(matchCluster).at(ii)->col < (*it).first )
           {
              posColToT.insert(it , std::make_pair(event.clusters.at(matchCluster).at(ii)->col, event.clusters.at(matchCluster).at(ii)->tot) );
	      found = true;
	      break;
           }
	   else if( event.clusters.at(matchCluster).at(ii)->col == (*it).first )
	   {
	      (*it).second+= event.clusters.at(matchCluster).at(ii)->tot;
	      //if(event.clusters.at(matchCluster).size() ==8 ) 
	      //{
	        clusterCol= event.clusters.at(matchCluster).at(ii)->col;
		//sharingToT= (*it).second;
	      //}
	      found = true;
	      pureRow = false;
	      break;
	   }
	   if(event.clusters.at(matchCluster).at(ii)->row != roww) pureRow=false;
	   //else
	   //{
	      //posColToT.push_back( std::make_pair(event.clusters.at(matchCluster).at(ii)->col, event.clusters.at(matchCluster).at(ii)->tot) );
	      //break;
	   //}
	}
	if(!found) posColToT.push_back( std::make_pair(event.clusters.at(matchCluster).at(ii)->col, event.clusters.at(matchCluster).at(ii)->tot) );
	
	found=false;
        for(std::vector<std::pair<int,int> >::iterator it=posRowToT.begin(); it!=posRowToT.end(); ++it)
	{
           if( event.clusters.at(matchCluster).at(ii)->row < (*it).first )
           {
              posRowToT.insert(it , std::make_pair(event.clusters.at(matchCluster).at(ii)->row, event.clusters.at(matchCluster).at(ii)->tot) );
	      found = true;
	      break;
           }
	   else if( event.clusters.at(matchCluster).at(ii)->row == (*it).first )
	   {
	      (*it).second+= event.clusters.at(matchCluster).at(ii)->tot;
	      found = true;
	      pureCol = false;
	      break;
	   }
	   if(event.clusters.at(matchCluster).at(ii)->col != colw) pureCol=false;
	   //if(posColToT.size()==7) pureCol=true;
	   //else
	   //{
	      //posRowToT.push_back( std::make_pair(event.clusters.at(matchCluster).at(ii)->row, event.clusters.at(matchCluster).at(ii)->tot) );
	      //break;
	   //}
	}
	if(!found) posRowToT.push_back( std::make_pair(event.clusters.at(matchCluster).at(ii)->row, event.clusters.at(matchCluster).at(ii)->tot) );
     }
      
  }
  
   for(unsigned int ii=0; ii < posColToT.size(); ++ii)
   {
       if( hh_matchPixelTotCol.find(posColToT.size()) == hh_matchPixelTotCol.end() || hh_matchPixelTotCol[posColToT.size()].find(ii) == hh_matchPixelTotCol[posColToT.size()].end() )
       {
          std::stringstream title;
          title << "Matching pixel charge [ToT] for pixel_" << ii;
          hh_matchPixelTotCol[posColToT.size()][ii] = new TH1D("", title.str().c_str(), 250, 0, 250);
       }
       if(hh_pureMatchPixelTotCol.find(posColToT.size()) == hh_pureMatchPixelTotCol.end() || hh_pureMatchPixelTotCol[posColToT.size()].find(ii) == hh_pureMatchPixelTotCol[posColToT.size()].end() )
	{
	    std::stringstream title; title << "Pure matching pixel charge [ToT] for pixel_" << ii;
            hh_pureMatchPixelTotCol[posColToT.size()][ii] = new TH1D("", title.str().c_str(), 16, 0, 16);
	}
	if(hh_shareMatchPixelTotCol.find(posColToT.size()) == hh_shareMatchPixelTotCol.end() || hh_shareMatchPixelTotCol[posColToT.size()].find(ii) == hh_shareMatchPixelTotCol[posColToT.size()].end() )
	{
	    std::stringstream title; title << "Share matching pixel charge [ToT] for pixel_" << ii;
            hh_shareMatchPixelTotCol[posColToT.size()][ii] = new TH1D("", title.str().c_str(), 50, 0, 50);
	}
	
       if(posRowToT.size()==2 && clusterCol==posColToT[ii].first)
       {
         hh_shareMatchPixelTotCol[posColToT.size()][ii]->Fill( posColToT[ii].second );
       }
       hh_matchPixelTotCol[posColToT.size()][ii]->Fill( posColToT[ii].second );
       if(pureRow) {hh_pureMatchPixelTotCol[posColToT.size()][ii]->Fill( posColToT[ii].second );}
   }   
   
   for(unsigned int ii=0; ii < posRowToT.size(); ++ii)
   {
       if( hh_matchPixelTotRow.find(posRowToT.size()) == hh_matchPixelTotRow.end() || hh_matchPixelTotRow[posRowToT.size()].find(ii) == hh_matchPixelTotRow[posRowToT.size()].end() )
       {
          std::stringstream title;
          title << "Matching pixel charge [ToT] for pixel in Row_" << ii;
          hh_matchPixelTotRow[posRowToT.size()][ii] = new TH1D("", title.str().c_str(), 250, 0, 250);
       }
       if(hh_pureMatchPixelTotRow.find(posRowToT.size()) == hh_pureMatchPixelTotRow.end() || hh_pureMatchPixelTotRow[posRowToT.size()].find(ii) == hh_pureMatchPixelTotRow[posRowToT.size()].end() )
	{
	    std::stringstream title; title << "Pure matching pixel charge [ToT] for pixel in Row_" << ii;
            hh_pureMatchPixelTotRow[posRowToT.size()][ii] = new TH1D("", title.str().c_str(), 250, 0, 250);
	}
       hh_matchPixelTotRow[posRowToT.size()][ii]->Fill( posRowToT[ii].second );
       if(pureCol) {hh_pureMatchPixelTotRow[posRowToT.size()][ii]->Fill( posRowToT[ii].second );}
   }
  
  if( doCharge == true ) { h_matchClusterQ->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster), event )); }
  
  h_matchClusterTotMap->Fill( cluster::getChargeWeightedCol(event.clusters.at(matchCluster)), cluster::getChargeWeightedRow(event.clusters.at(matchCluster)), cluster::getSumTot(event.clusters.at(matchCluster)) );
  if( doCharge == true ) { h_matchClusterQMap->Fill( cluster::getChargeWeightedCol(event.clusters.at(matchCluster)), cluster::getChargeWeightedRow(event.clusters.at(matchCluster)), cluster::getSumChargePP(event.clusters.at(matchCluster),event) ); }
  
  // Not sure what this is checking for
  // and cut doesn't generalize to FE-I4 (BJD)
  if( tbutils::getFoldedX(event) < 300 ) {
    if( doCharge == true ) { h_clusterq_no_bg->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster),event )); }
  }

  switch( event.clusters.at(matchCluster).size() ) {
  case 1:
    if( doCharge == true ) { h_matchClusterQ_1->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster),event )); }
    break;
  case 2:
    if( doCharge == true ) { h_matchClusterQ_2->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster),event )); }
    break;
  case 3:
    if( doCharge == true ) { h_matchClusterQ_3->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster),event)); }
    break;
  case 4:
    if( doCharge == true ) { h_matchClusterQ_4->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster),event)); }
    break;
  default:
    if( doCharge == true ) { h_matchClusterQ_5plus->Fill( cluster::getSumChargePP( event.clusters.at(matchCluster),event )); }
    break;
  }
  
//  h_mapClusterChargePP->Fill(cluster::getEtaCorrectedCol( event.clusters.at(matchCluster) event.dut ),
//			cluster::getEtaCorrectedRow( event.clusters.at(matchCluster), event.dut ),
//			cluster::getSumChargePP( event.clusters.at(matchCluster),event ) ); 

  // Use only 1-hit clusters
  //if( event.clusters.at(matchCluster).size() != 1 ) {
  //  TBALOG(kDEBUG2) << "Use only 1-hit clusters for ToT in electrode, skip this event with cluster size = " << event.clusters.at(matchCluster).size() << endl;
   // return;
 // }
  
  // Find the distance to the electrode
  double elecEdgeDistX = tbutils::getTrackElectrodeDistance(event.trackX, event.dut->getPitchX(), event.dut->getNumElec());
  double elecEdgeDistY =  tbutils::getTrackElectrodeDistance(event.trackY, event.dut->getPitchY(), 1);

  // Move the coordinates so that they correspond to the expected center of the electrode
  TVector2 eVec(elecEdgeDistX - 0.5*(event.dut->getPitchX()/event.dut->getNumElec()), elecEdgeDistY - 0.5*(event.dut->getPitchY()/1.0));
  h_elecMap->Fill(eVec.X(), eVec.Y());
  h_elecMapTot->Fill(eVec.X(), eVec.Y(), cluster::getSumTot(event.clusters.at(matchCluster)));
  h_elecDist->Fill(eVec.Mod());
  h_elecDistTot->Fill(eVec.Mod(), cluster::getSumTot(event.clusters.at(matchCluster)));
  
}


void SumTot::finalize(const TbConfig &config){

  config.drawAndSave(name, "fullChipTot", h_fullChipTot);
  config.drawAndSave(name, "fullChipQ", h_fullChipQ);
  config.drawAndSave(name, "maxCellTot", h_maxCellTot);
  config.drawAndSave(name, "maxCellQ", h_maxCellQ);
  config.drawAndSave(name, "maxClusterTot", h_maxClusterTot);
  config.drawAndSave(name, "maxClusterQ", h_maxClusterQ);
  config.drawAndSave(name, "matchClusterTot", h_matchClusterTot);
  config.drawAndSave(name, "matchPixelTot", h_matchPixelTot);
  
  for(std::map< int, std::map<int,TH1D*> >::iterator siz=hh_matchPixelTotCol.begin(); siz!=hh_matchPixelTotCol.end(); ++siz)
  {
     for(std::map<int,TH1D*>::iterator it=siz->second.begin(); it!=siz->second.end(); ++it)
     { 
        std::stringstream title;
        title << "matchPixelTotCol_" << siz->first << "-" << it->first;
        config.drawAndSave(name, title.str().c_str(), it->second);
	exit(12);
     }
  }
  for(std::map< int, std::map<int,TH1D*> >::iterator siz=hh_pureMatchPixelTotCol.begin(); siz!=hh_pureMatchPixelTotCol.end(); ++siz)
  {
     for(std::map<int,TH1D*>::iterator it=siz->second.begin(); it!=siz->second.end(); ++it)
     { 
        std::stringstream title;
        title << "pureMatchPixelTotCol_" << siz->first << "-" << it->first;
        config.drawAndSave(name, title.str().c_str(), it->second);
     }
  }  
  for(std::map< int, std::map<int,TH1D*> >::iterator siz=hh_shareMatchPixelTotCol.begin(); siz!=hh_shareMatchPixelTotCol.end(); ++siz)
  {
     for(std::map<int,TH1D*>::iterator it=siz->second.begin(); it!=siz->second.end(); ++it)
     { 
        std::stringstream title;
        title << "shareMatchPixelTotCol_" << siz->first << "-" << it->first;
        config.drawAndSave(name, title.str().c_str(), it->second);
     }
  }
  
  for(std::map< int, std::map<int,TH1D*> >::iterator siz=hh_matchPixelTotRow.begin(); siz!=hh_matchPixelTotRow.end(); ++siz)
  {
     for(std::map<int,TH1D*>::iterator it=siz->second.begin(); it!=siz->second.end(); ++it)
     { 
        std::stringstream title;
        title << "matchPixelTotRow_" << siz->first << "-" << it->first;
        config.drawAndSave(name, title.str().c_str(), it->second);
     }
  }
  for(std::map< int, std::map<int,TH1D*> >::iterator siz=hh_pureMatchPixelTotRow.begin(); siz!=hh_pureMatchPixelTotRow.end(); ++siz)
  {
     for(std::map<int,TH1D*>::iterator it=siz->second.begin(); it!=siz->second.end(); ++it)
     { 
        std::stringstream title;
        title << "pureMatchPixelTotRow_" << siz->first << "-" << it->first;
        config.drawAndSave(name, title.str().c_str(), it->second);
     }
  }
  config.drawAndSave(name, "matchClusterQ", h_matchClusterQ);
  config.drawAndSave(name, "matchClusterQ_size1", h_matchClusterQ_1);
  config.drawAndSave(name, "matchClusterQ_size2", h_matchClusterQ_2);
  config.drawAndSave(name, "matchClusterQ_size3", h_matchClusterQ_3);
  config.drawAndSave(name, "matchClusterQ_size4", h_matchClusterQ_4);
  config.drawAndSave(name, "matchClusterQ_size5plus", h_matchClusterQ_5plus);
  config.saveToFile(name, "clusterq_noBiasGrid", h_clusterq_no_bg); // BJD: ?

  config.drawAndSave(name, "elecMap", h_elecMap);
  config.drawAndSave(name, "elecMapToT", h_elecMapTot);
  config.drawAndSave(name, "elecDist", h_elecDist);
  config.drawAndSave(name, "elecDistTot", "col", h_elecDistTot);

  config.drawAndSave(name, "matchClusterTotMap", h_matchClusterTotMap);
  totMap(config, name, "totMap", h_matchClusterTotMap);
  
  totElecMap(config, name, "totElecMap", h_elecMapTot);  
  ElecDistTot(config, name, "elecDistTot", h_elecDistTot );
  
}


void SumTot::totMap( const TbConfig& config, const char* analysisName, const char* histoName, TH3D* totmap3d ) const {

  char* fileName = config.buildHistName(analysisName,histoName);

  // plot the mean tot in each pixel
  TProfile2D * totmap2d = totmap3d->Project3DProfile("yx");
  totmap2d->SetTitle("Matching Cluster ToT;Column;Row");
  //totmap2d->SetName(TString::Format("totmap2d_%s",histoName));

  TCanvas* canvas = new TCanvas();
  canvas->cd();
  totmap2d->Draw("colz");
  canvas->SaveAs(fileName);  

  delete fileName;
  delete totmap2d;
  delete canvas;
  
}


void SumTot::totElecMap( const TbConfig& config, const char* analysisName, const char* histoName, TH3D* totmap3d ) const {

  char* fileName = config.buildHistName(analysisName,histoName);

  // Plot the mean tot in each pixel
  TProfile2D * totmap2d =  totmap3d->Project3DProfile("yx");
  totmap2d->SetStats(false);
  totmap2d->SetTitle(";Distance from read-out electrode center [#mum];Distance from read-out electrode center [#mum];Matching Cluster ToT");
  //totmap2d->SetName(TString::Format("totmap2d_%s",histoName));

  TCanvas* canvas = new TCanvas();
  canvas->cd();
  totmap2d->Draw("cont4z");
  canvas->SaveAs(fileName);

  delete fileName;
  delete totmap2d;
  delete canvas;
  
}


void SumTot::ElecDistTot( const TbConfig& config, const char* analysisName, const char* histoName, TH2D* elecDistTot ) const {
  
  //project each bin
  elecDistTot->Rebin2D(2,1);
  int nBins = elecDistTot->GetNbinsX();  

  // Define two distances that should be compared directly
  double d1 = 20.0;  // outside electrode!
  double d2 = 4.0; // inside electrode!
  TString d1s = "";
  TString d2s = "";
  TH1D* h1 = 0;
  TH1D* h2 = 0;
  
  for( int i = 1; i != nBins; ++i) {
    TH1D* h = (TH1D*) elecDistTot->ProjectionY(TString::Format("pr%i",i),i,i,"e");
    h->SetStats(false);
    h->SetTitle(";Cluster charge [ToT];Tracks");

    TCanvas* canvas = new TCanvas();
    canvas->cd();
    //h->SetFillColor(kOrange);
    h->SetFillColor(kRed-4);
    TH1D* h_copy = (TH1D*) h->DrawCopy("hist");
    double xMin(0), xMax(0);
    xMin = h_copy->GetXaxis()->GetBinLowEdge(i);
    xMax = h_copy->GetXaxis()->GetBinUpEdge(i);
    tbutils::myText(0.15, 0.85, kBlack, TString::Format("[%.1f,%.1f] #mum",xMin,xMax));
    
    TString fname = TString::Format("%s%i",histoName,i);
    char* fileName = config.buildHistName(analysisName,fname.Data());
    canvas->SaveAs(fileName);

    config.saveToFile(this->name, fname, h);  
    delete fileName;
    delete canvas;
    //delete h_copy;
    
    if( d1>xMin and d1<=xMax) {
      h1 = h;
      d1s = TString::Format("[%.1f,%.1f] #mum",xMin,xMax);
    } else if( d2>xMin and d2<=xMax ) {
      h2 = h;
      d2s = TString::Format("[%.1f,%.1f] #mum",xMin,xMax);
    }
  }

  if( h1 and h2 ) {
    TCanvas* canvas_Comp = new TCanvas();
    canvas_Comp->cd();   
    h2->SetTitle(";Arbitrary Units;Cluster charge [ToT]");
    h2->SetFillColor(kRed-4);
    h2->SetFillStyle(1001);      
    h2->DrawCopy("hist");
    //tbutils::Norm(h1, 1.0, 1); // BJD: Something was going wrong here
    //tbutils::Norm(h2, 1.0, 1); // so I commented it out :)
    h1->SetFillColor(kBlack);
    h1->SetFillStyle(3004);      
    h1->DrawCopy("hist same");

    TLegend* leg = new TLegend(0.65,0.7,0.88,0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0); 
    leg->AddEntry(h1, d1s, "F");
    leg->AddEntry(h2, d2s, "F");
    leg->Draw();

    TString fname = TString::Format("%sComp",histoName);
    char* fileName = config.buildHistName(analysisName,fname.Data());
    canvas_Comp->SaveAs(fileName);  
    //canvas_Comp->SaveAs("temp.C");  
    delete fileName;
    delete canvas_Comp;
    delete leg;
  } 
  //delete elecDistTot;
 
}
