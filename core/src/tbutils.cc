#include "tbutils.h"
#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TObjString.h"
#include "TPave.h"
#include <assert.h>

using std::cout;
using std::endl;

bool tbutils::inCentralRegion(PllHit* hit, const Event& event){
  bool accept(true);
  if( hit->col < event.dut->getSkipCols() or 
      hit->col > (event.dut->getMaxCol() - event.dut->getSkipCols())){ accept = false;}
  if( hit->row < event.dut->getSkipRows() or 
      hit->row > (event.dut->getMaxRow() - event.dut->getSkipRows())){ accept = false;}
  return(accept);
}

// bool tbutils::trackInCentralRegion(const double& trackX,const double& trackY) {
//   double X = trackX + 0.5*event::pitchX;
//   double Y = trackY + 0.5*event::pitchY;      
//   if( X < (event::pitchX * event::skipCols) ) return false;
//   if( X > event::pitchX * (event::maxCol - event::skipCols) ) return false;
//   if( Y < (event::pitchY * event::skipRows) ) return false;
//   if( Y > event::pitchY * (event::maxRow - event::skipRows) ) return false;
//   return true;
// } 

bool tbutils::isMatch(const PllHit* const hit, const Event &event){
  double residualX = fabs( event.trackX - (event.dut->getPitchX() * hit->col) );
  double residualY = fabs( event.trackY - (event.dut->getPitchY() * hit->row) );
  if( residualX < event.dut->getMatchX() and 
      residualY < event.dut->getMatchY() ) {
    return(true);
  } else {
    return(false);
  }
}

int tbutils::getCol(double trackX, const Event &event){
  int col = (int) ((trackX + 0.5 * event.dut->getPitchX()) / event.dut->getPitchX());
  if(col < 0 or col > event.dut->getMaxCol()){ col = -1; }
  return( col );
}

int tbutils::getRow(double trackY, const Event &event){
  int row = (int) ((trackY + 0.5 * event.dut->getPitchY()) / event.dut->getPitchY());
  if( row < 0 or row > event.dut->getMaxRow()){ row = -1; }
  return( row );
}

double tbutils::getTrackEdgeDistance(const double& trackX, const double& p ) {
  //Find the distance of the track to the closest edge of the pixel 
  //shift the track coordinate origo to (0,0) to be able to use modulo
  double trackXT = trackX + p*0.5;
  while(trackXT > p){ trackXT -= p; }
  if(trackXT > 0.5*p) trackXT = p - trackXT;
  return trackXT;
}
double tbutils::getTrackEdgeDistance(const Event &event, TString side) {
  side.ToLower();
  return side == "x" ? getTrackEdgeDistance(event.trackX,event.dut->getPitchX()) : getTrackEdgeDistance(event.trackY,event.dut->getPitchY());
}
double tbutils::getTrackTwoPixelEdgeDistance(const double& trackX, const double& p ) {
  //gets the distance from a point 1/2 pixel off and extending over 2 pixel dimension
  double trackXT = trackX;
  while(trackXT > 2*p){ trackXT -= 2*p; }
  return trackXT;
}
double tbutils::getTrackTwoPixelEdgeDistance(const Event &event, TString side) {
  side.ToLower();
  return side == "x" ? getTrackTwoPixelEdgeDistance(event.trackX,event.dut->getPitchX()) : getTrackTwoPixelEdgeDistance(event.trackY,event.dut->getPitchY());
}

double tbutils::getTrackElectrodeDistance(const double& trackX, const double& p, const int& nelectrodes) {
  //Find the track coordinates wrt to the area around the electrode 
  //First get the distance to the edge
  double trackXT = trackX + p*0.5;
  while(trackXT > p){ trackXT -= p; }
  //To get the box centered on the electrode we need to subtract
  //the electrode distance 
  double electrode_spacing = p/static_cast<double>(nelectrodes);
  while(trackXT > electrode_spacing){ trackXT -= electrode_spacing; }
  return trackXT;
}

double tbutils::getFoldedX(const Event& event) {
	double x = fmod(event.trackX+0.5*event.dut->getPitchX(), 2*event.dut->getPitchX());
	if (x > event.dut->getPitchX()) x = 2*event.dut->getPitchX() - x;
	return x;
}

double tbutils::getFoldedY(const Event& event) {
	double y = fmod(event.trackY+0.5*event.dut->getPitchY(), 2*event.dut->getPitchY());
	if (y > event.dut->getPitchY()) y = 2*event.dut->getPitchY() - y;
	return y;
}

// added by Botho 2014-10-20
double tbutils::getPixelX(const Event& event) {
	return fmod(event.trackX+0.5*event.dut->getPitchX(), event.dut->getPitchX());
}

double tbutils::getPixelY(const Event& event) {
	return fmod(event.trackY+0.5*event.dut->getPitchY(), event.dut->getPitchY());
}

void tbutils::dumpEventInfo(const Event& event) {
  if(event.track == NULL){
    cout << "The function dumpEventInfo is intended for BAT data, but cannot find Track object." << endl;
    return;
  }

  const Track& track = *event.track;
  cout << "\n----Track dump" << endl;
  cout << setw(20) << "trig " << track.trig << endl;
  cout << setw(20) << "chi2 " << track.ct << endl;
  cout << setw(20) << "chi2 (weighted)" << track.gt << endl;
  cout << setw(20) << "hitcount " << track.ft << endl;
  cout << "--Track params " << endl;
  for(int par = 0; par < track.nTrackParams; par++) { 
    TrackParams* params = (TrackParams*)track.trackParams[par]; 
#ifdef FOUR_PARAMS
    TString s  = TString::Format("%f,%f,%f,%f", 
				 params->params[0] * 1000.0,
				 params->params[1] * 1000.0,
				 params->params[2],
				 params->params[3]);
#else
    TString s  = TString::Format("%f,%f,%f,%f,%f", 
				 params->params[0] * 1000.0,
				 params->params[1] * 1000.0,
				 params->params[2],
				 params->params[3],
				 params->params[4]);
#endif
    cout << "iden " << params->iden << ": " << s.Data() << endl;
  }
  cout << "\n----PllHit dump" << endl;
  cout << "hits " << track.nPllHit << endl;
  cout << "(iden,row,col,tot)" << endl;  
  for(int ii = 0; ii < track.nPllHit; ++ii){
    PllHit* hit = (PllHit*)track.pllHit[ii];
    if(!hit) continue;
    cout << hit->iden << ","<<hit->row << ","<<hit->col << ","<<hit->tot << endl;
  }
}

bool tbutils::NAND(bool a, bool b) {
  return (!a||!b);
}

bool tbutils::XOR(bool a, bool b) {
  bool n = NAND(a,b);
  return NAND(NAND(a,n),NAND(b,n));
}


double tbutils::turnon_v2(Double_t *x, Double_t *par)
{
  double mu = par[0];
  double sigma = par[1];
  double plateau = par[2];
  //double offset = par[3];
  double offset = 0;
  double x0 = x[0];
  
  double arg = 0;
  arg = (x0 - mu)/(sqrt(2.0)*sigma);
  //arg = (pt - halfpoint)/ (TMath::Sqrt(pt)*slope);
  double fitval = offset+0.5*plateau*(1+TMath::Erf(arg));
  return fitval;

}



//TF1* HitEff::GetGausLineFit(TH1D* h,TString s,double xmin,double xmax,double hp,double sl,double pla) {
TF1* tbutils::GetGausLineFit(TObject* h,TString s,double xmin,double xmax,
			     double hp,double sl,double pla) {
  TF1 *f = new TF1(TString::Format("turnon_%s",s.Data()),turnon_v2,xmin,xmax,3);
  f->SetParName(0,"halfpoint");
  f->SetParName(1,"slope");
  f->SetParName(2,"plateau");
  f->SetParameter(0,hp);
  f->SetParameter(1,sl);
  f->SetParameter(2,pla);
  TGraph* gr=0;
  TH1D* hh=0;
  if(h->InheritsFrom("TGraph")) gr=(TGraph*)h;
  if(h->InheritsFrom("TH1D")) hh = (TH1D*)h;
  assert(XOR(gr,hh));
  if(hh) hh->Fit(f,"0R"); //  h->Fit(f,"V0R");
  if(gr) gr->Fit(f,"0R"); //  h->Fit(f,"V0R");
  return f;
}

void tbutils::myText(Double_t x,Double_t y,Color_t color,const char *text) {  
  //Double_t tsize=0.05;
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

void tbutils::myMultiLineText(Double_t x1,Double_t y1, TString text,
			      Color_t boxcolor, Color_t textcolor) {
  
  Double_t textsizeX=0.017;
  Double_t textsizeY=0.044;
  
  //split up in lines
  TObjArray* lines = text.Tokenize('|');      
  Int_t nlines=lines->GetEntries();
  int nxmax=0;
  //std::cout << "nlines="<<nlines<<" nxmax="<<nxmax<<std::endl;
  for(int i=0;i!=nlines;++i) {            
    const TObjString* lineobj = (TObjString*)lines->At(i);
    if(!lineobj) {std::cout<<"no lineobj!"<<std::endl; exit(0);}
    TString line=lineobj->GetString();
    //std::cout << "line "<<i<<"="<<line.Data()<<std::endl;
    Int_t ll=line.Length();
      if(ll>nxmax) nxmax=ll;
  }
  //std::cout << "nlines="<<nlines<<" nxmax="<<nxmax<<std::endl;
  Double_t y1new=y1-0.025;
  Double_t x1new=x1-0.01;
  Double_t y2=y1-nlines*textsizeY;
  Double_t x2=x1+nxmax*textsizeX;
  
  TPave mbox(x1new,y1new,x2,y2);
  //mbox->SetShadowColor(0);
   //mbox->SetFillColor(boxcolor);
   //mbox->SetFillStyle(1001);
  mbox.Draw();
  
  for(int i=0;i!=nlines;++i) {
    const TObjString* lineobj = (TObjString*)lines->At(i);
    if(!lineobj) {std::cout<<"no lineobj!"<<std::endl; exit(0);}
    TString line=lineobj->GetString();      
    //std::cout << "i="<<i<<" line="<<line.Data()<<std::endl;
    TLatex l;
    l.SetTextColor(textcolor);
    l.SetTextAlign(12); //l.SetTextSize(tsize); 
    l.SetNDC();          
    l.DrawLatex(x1,y1-double(i)*textsizeY,line.Data());          
      //std::cout << "drawing line="<<line.Data()<<std::endl;
  }
  if(lines) {
    delete lines;
  }
  return;  
}


void tbutils::atlasHistStyle() {  
  //std::cout << std::endl << "Welcome to the ATLAS rootlogon.C" << std::endl;
  //
  // based on a style file from BaBar
  //
  //..BABAR style from RooLogon.C in workdir
  TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");
  
  // use plain black on white colors
  Int_t icol=0;
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetStatColor(icol);
  atlasStyle->SetFillColor(icol);
  
  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);
  atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.12);
  
  // use large fonts
  //Int_t font=72;
  Int_t font=42;
  Double_t tsize=0.05;
  atlasStyle->SetTextFont(font);
  
  
  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");
  
  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");
  
  
  //use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth((Width_t)2.);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  
  //get rid of X error bars and y error bar caps
  //atlasStyle->SetErrorX(0.001);
  
  //do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);
  
  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);
  
  //gROOT->SetStyle("Plain");
  
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  
  //std::cout << "... finished!...type: $>\'gROOT->SetStyle(\"ATLAS\"); \'gROOT->ForceStyle();\'" << std::endl;
  gROOT->SetStyle("ATLAS"); 
  gROOT->ForceStyle();
   
}

void tbutils::HistStyle() {  
  //std::cout << std::endl << "Welcome to the ATLAS rootlogon.C" << std::endl;
  //
  // based on a style file from BaBar
  //
  //..BABAR style from RooLogon.C in workdir
  //TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");
  
  // use plain black on white colors
  Int_t icol=0;
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetStatColor(icol);
  gStyle->SetFillColor(icol);
  
  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.12);
  
  // use large fonts
  //Int_t font=72;
  Int_t font=42;
  Double_t tsize=0.05;
  gStyle->SetTextFont(font);
  
  
  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");
  
  gStyle->SetLabelSize(tsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(tsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(tsize,"z");
  gStyle->SetTitleSize(tsize,"z");
  
  
  //use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth((Width_t)2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  
  //get rid of X error bars and y error bar caps
  //gStyle->SetErrorX(0.001);
  
  //do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);
  
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  //gROOT->SetStyle("Plain");
  
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  
  //std::cout << "... finished!...type: $>\'gROOT->SetStyle(\"ATLAS\"); \'gROOT->ForceStyle();\'" << std::endl;
  //gROOT->SetStyle("ATLAS"); 
  gROOT->ForceStyle();
   
}

double tbutils::fitFunc(double *v, double *par){
  //Two error functions to fit a box 
  double x=v[0];
  double sigma=par[0];
  double area=par[1];
  double x0=par[2];
  double width=par[3];
  double angCoeff=par[4];

  double sfita = 0.5*(1+TMath::Erf((x+width-x0)/(TMath::Sqrt(2.)*sigma)));
  double sfitb = 0.5*(1+TMath::Erf((x-width-x0)/(TMath::Sqrt(2.)*sigma)));
  
  return 0.5*(1 + angCoeff*(x-x0))*area*(sfita -sfitb)/(width);
}

double tbutils::fitFuncBox(double *v, double *par){
  //Two error functions to fit a box 
  double x=v[0];
  double sigma=par[0];
  double area=par[1];
  double width=par[2];

  double sfita = 0.5*(1+TMath::Erf((x+width)/(TMath::Sqrt(2.)*sigma)));
  double sfitb = 0.5*(1+TMath::Erf((x-width)/(TMath::Sqrt(2.)*sigma)));
  
  return 0.5*(1 + area*(sfita -sfitb)/(width));
}

double tbutils::multiG(double *v, double *par){

  double x = v[0];

  double coreMean = par[0];
  double outMean = par[1];

  double coreWidth = par[2];
  double outWidth = par[3];

  double coreFrac = par[4];

  double norm = par[5];

  // G = norm * [coreFrac * C + (1 - coreFrac) * O ]

  // C: core gaussian
  double coreArg = (x-coreMean);
  coreArg /= coreWidth;
  coreArg *= coreArg;
  coreWidth /= 2.0;

  double C = TMath::Exp(-coreArg);
  C /= (TMath::Sqrt(2.0*TMath::Pi())*coreWidth);


  // O: outlier gaussian
  double outArg = (x-outMean);
  outArg /= outWidth;
  outArg *= outArg;
  outWidth /= 2.0;

  double O = TMath::Exp(-outArg);
  O /= (TMath::Sqrt(2.0*TMath::Pi())*outWidth);



  // G = norm *[ coreFrac * C + (1 - coreFrac) * O ]
  double mG = norm*(coreFrac * C + (1. - coreFrac)*  O );

  return mG;
}

void tbutils::Norm(TH1* histo, double norm, Int_t IncludeOUBin){
  if(IncludeOUBin) {//default value      
    if (histo->Integral(0,histo->GetNbinsX()+1))
      histo->Scale( norm / histo->Integral(0,histo->GetNbinsX()+1));
    else
      std::cout << "\n\nWARNING: integral is zero...normalization does not work!\n" << std::endl;
  }
  else{
    if (histo->Integral(1,histo->GetNbinsX()))
      histo->Scale( norm / histo->Integral(1,histo->GetNbinsX()));
    else
      std::cout << "\n\nWARNING: integral is zero...normalization does not work!\n" << std::endl;
  }
}
