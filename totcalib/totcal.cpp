//#include "/afs/cern.ch/user/p/phansson/scripts/root_macros/gen_histos_light.hpp"
#include "rootlogonCompile.h"
//#include "utilsCompile.h"
#include "TGraphErrors.h"

#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TObjString.h"
#include "TPave.h"

using std::cout;
using std::endl;

double Cinj;
TFile* outfile;
std::string name;
int debug = 0;
std::list<int*> masks;

void myBoxMultiLineText(Double_t x1, 
                        Double_t y1,
                        TString text,
                        Color_t boxcolor=kWhite,
                        Color_t textcolor=kBlack) {
   
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

void myText(Double_t x,Double_t y,Color_t color,const char *text) {

//   //Double_t tsize=0.05;
   TLatex l; l.SetTextAlign(12);// l.SetTextSize(tsize); 
   l.SetNDC();
   l.SetTextColor(color);
   l.DrawLatex(x,y,text);
}


double totqfnc(double* x, double* par) {
  //(TOT = p0 (p1+Q)/(p2+Q))
  double Q = x[0];
  double tot = par[0]*(par[1]+Q)/(par[2]+Q);
  return tot;
}

bool IsMasked(const int& row, const int& col) {
  for(std::list<int*>::iterator it = masks.begin(); it != masks.end(); it++){
    if( (col == (*it)[0]) && (row == (*it)[1])){ return true;}
  }
  return false;
}

void addMask(const int& col, const int& row){
  int* mask = new int[2];
  mask[0] = col; mask[1] = row;
  masks.push_back(mask);
}

void addMask(const std::string&  file) {
  /* Use a text file with a list of row,column as a mask
   */
  if(file=="") {
    cout << "Warning: no mask file exist for module " << name << " !!!" << endl;
    return;
  }
  //the file should contain a list of pixels with 'row,column' with line breaks
  std::string line;
  std::ifstream fs(file.c_str());
  if(!fs.is_open()) {
    cout << "The mask '" << file << "' could not be opened!!" << endl;
    exit(0);         
  } else {
    while(!fs.eof()) {
      getline(fs,line);
      std::string::size_type found=line.find(":");
      if (found==std::string::npos) {
	//cout << "This mask line have wrong configuration: " << line << endl;
	//exit(0);         
	continue;
      }
      std::string scol = line.substr(found+1);
      std::string srow = line.erase(found,std::string::npos);
      std::istringstream issc(scol);
      int column;
      issc >> column;
      std::istringstream issr(srow);
      int row;
      issr >> row;
      addMask(column,row);
    }
  }  
  fs.close();
  return;
}

/*! \brief Helper class tot cal --> root conversion standalone program.
 *
 */
struct pixel {
  int col;
  int row;
  std::map<int,TH1F*> histos;
  std::map<int,TF1*> tot_fit;
  TGraphErrors* tot_cal;
  TF1* tot_cal_fit;
  pixel(const int& r,const int& c) { col=c; row=r; }
  void Fill(const int& vcal,const int& tot, const int& n=1) {
    std::map<int,TH1F*>::iterator ih = histos.find(vcal);
    if(ih!=histos.end()) {
      if(debug) cout<<"filling vcal="<<vcal<<" tot="<<tot<<" n="<<n<<std::endl;
      for(int i=0;i!=n;i++) (ih->second)->Fill(tot);
    } else {
      //std::cout<<"ERROR! this "<<vcal<<" hist is not created!"<<std::endl;exit(0);
      float q = float(Cinj)*TMath::Power(10,-15)*float(vcal)*TMath::Power(10,-3)/(1.6*TMath::Power(10,-19));
      float k = 60.0/20000.0;
      int appr_tot = TMath::Nint(q*k);//TMath::FloorNint(q*k);
      if(name=="160") appr_tot *= 0.8;
      float xmin = appr_tot-20;
      float xmax = appr_tot+20;
      int N = 40;
      TH1F* h = new TH1F(TString::Format("h_tot_%i_%i_%i",row,col,vcal),TString::Format("h_tot_%i_%i_%i;ToT;Entries",row,col,vcal),N,xmin,xmax);
      if(debug) cout<<"vcal="<<vcal<<" q="<<q<<" a_tot="<<appr_tot<<" xmin="<<xmin<<" xmax="<<xmax<<endl;
      histos[vcal] = h;
      this->Fill(vcal,tot,n);
    }
  }
  
  void FitTot(const int& vcal) {    
    std::map<int,TH1F*>::iterator ih = histos.find(vcal);
    if(ih!=histos.end()) {
      (ih->second)->Fit("gaus","Q");  
      TF1* fnc = (ih->second)->GetFunction("gaus");
      tot_fit[vcal] = (TF1*)fnc->Clone(TString::Format("f_tot_fit_%i_%i_%i",row,col,vcal));
    }
  }
  
  
  void MakeCal() {
    TGraphErrors* gr = new TGraphErrors((int)tot_fit.size()); 
    double m,em;
    int dummy=0;
    for(std::map<int,TF1*>::iterator i=tot_fit.begin();i!=tot_fit.end();++i) {
      m = i->second->GetParameter(1);
      em = i->second->GetParError(1);
      float q = float(Cinj)*TMath::Power(10,-15)*float(i->first)*TMath::Power(10,-3)/(1.6*TMath::Power(10,-19));
      gr->SetPoint(dummy,q,m);
      gr->SetPointError(dummy,q*0.01,em);
      dummy++;
    }
    
    TF1* f = new TF1(TString::Format("fcal_%i_%i",row,col),totqfnc,4000,45000,3);    
    f->SetParameter(0,1000);
    f->SetParameter(1,-2000);
    f->SetParameter(2,200000);    
    f->SetTitle(TString::Format("ToT Calbration;Injected charge [e^{-}];ToT"));
    gr->Fit(f,"QR");
    tot_cal=gr;
    tot_cal_fit=f;
  }
  

  void PlotTot(const int& vcal) {
    TCanvas* c = new TCanvas("c","c",10,20,700,500);
    std::map<int,TH1F*>::iterator ih = histos.find(vcal);
    if(ih!=histos.end()) {
      (ih->second)->Draw();  
      c->SaveAs(TString::Format("tot_%i_%i_%i_%s.png",row,col,vcal,name.c_str()));
      delete c;
    }
  }
  void PlotCal() {
    TCanvas* c = new TCanvas("c","c",10,20,700,500);
    tot_cal->Draw("AP");  
    myText(0.75,0.27,kBlack,TString::Format("row %i col %i",row,col));
    myText(0.75,0.22,kBlack,TString::Format("%s",name.c_str()));
    myBoxMultiLineText(0.15,0.83,TString::Format("p_{0}=%.1f|p_{1}=%.1f|p_{2}=%.1f",
						 tot_cal_fit->GetParameter(0),
						 tot_cal_fit->GetParameter(1),
						 tot_cal_fit->GetParameter(2)));    
    c->SaveAs(TString::Format("totcal_%i_%i_%s.png",row,col,name.c_str()));
    delete c;
  }
};

/*! \brief Helper class tot cal --> root conversion standalone program.
 *
 */
struct pixels {
  std::vector<pixel> v;
  pixel* FindPixel(const int& row,const int& col) {
    std::vector<pixel>::iterator i=v.begin();
    std::vector<pixel>::iterator iE=v.end();
    for(;i!=iE;++i) {if(i->row==row && i->col==col) return &(*i);}
    return NULL;
  }
  pixel* Add(const int& row,const int& col) {
    pixel p(row,col);
    v.push_back(p);
    return &(v.back());
  }

  void SaveAll(TFile* file) {
    TDirectory* dir = gDirectory;
    file->cd();
    for(std::vector<pixel>::iterator i=v.begin();i!=v.end();++i) {      
      i->tot_cal_fit->Write(TString::Format("ftotToQ_%i_%i_%s",i->row,i->col,name.c_str()));
    }
    dir->cd();
  }

  void MakeCalAll() {
    int dummy=0;
    for(std::vector<pixel>::iterator i=v.begin();i!=v.end();++i) {      
      for(std::map<int,TH1F*>::iterator iv=(*i).histos.begin();iv!=(*i).histos.end();++iv) {
	int vcal = iv->first;
	(*i).FitTot(vcal);
      }//iv
      ++dummy;
      if(i->row==0) cout<<" Doing fit & cal for column "<<i->col<<endl;
      (*i).MakeCal();
    }//i
  }

  void PlotCalAll() {

    TCanvas* c = new TCanvas("c","c",10,20,700,500);      
    int dummy=0;
    for(std::vector<pixel>::iterator i=v.begin();i!=v.end();++i) {      
      TF1* f= (*i).tot_cal_fit;
      if(f) {
	f->SetLineColor(dummy+1);
	f->SetLineWidth(0.7);
	if(dummy==0) {
	  f->Draw();
	  TH1* h_0 = f->GetHistogram();	
	  h_0->SetMaximum(120.);
	}
	else {f->Draw("same");}
	dummy++;
	if(dummy%200==0)cout<<" plotted"<<dummy<<" cal fits"<<endl;
      }
    }
    c->SaveAs(TString::Format("totcalfitall_%s.png",name.c_str()));      
    delete c;
    

    std::vector<int> vvcal;
    for(std::vector<pixel>::iterator i=v.begin();i!=v.end();++i) {      
      for(std::map<int,TH1F*>::iterator iv=(*i).histos.begin();iv!=(*i).histos.end();++iv) {
	vvcal.push_back(iv->first);
      }
      if(vvcal.size()>0) break;
    }
    
    for(std::vector<int>::iterator iv=vvcal.begin();iv!=vvcal.end();++iv) {
      int vcal = *iv;
      c = new TCanvas("c","c",10,20,700,500);      
      dummy=0;
      for(std::vector<pixel>::iterator i=v.begin();i!=v.end();++i) {      
	TF1* f= (*i).tot_fit[vcal];
	if(f) {
	  f->SetLineColor(dummy+1);
	  f->SetLineWidth(0.7);
	  if(dummy==0) {
	    f->Draw();
	    TH1* h_0 = f->GetHistogram();	
	    h_0->SetMaximum(29.);	    
	  } 
	  else {f->Draw("same");}
	  dummy++;
	}
      }//i
      myText(0.75,0.27,kBlack,TString::Format("VCal=%imV",vcal));
      myText(0.75,0.22,kBlack,TString::Format("%s",name.c_str()));
      c->SaveAs(TString::Format("totcalall_%i_%s.png",vcal,name.c_str()));      
      delete c;
      cout<<" plotted ToT fit for all pixel for VCal="<<vcal<<endl;
      
    }//iv    
    
  }
  
};



int PixelType(int row,int col){  
  std::pair<int,int> pr = std::make_pair(row,col);
  if (!(pr.second%18 == 0 || pr.second%18 == 17)) {
    if (pr.first==152 || pr.first==154 || pr.first==156 || pr.first==158 || pr.first==161 || pr.first==163 || pr.first==165 || pr.first==167)
      return 4; // inter-ganged pixel
    else if (pr.first==153 || pr.first==155 || pr.first==157 || pr.first==159 || pr.first==160 || pr.first==162 || pr.first==164 || pr.first==166)
      return 2; // ganged pixel 
    else
      return 0; // normal pixel
  } else {
    if (pr.first==152 || pr.first==154 || pr.first==156 || pr.first==158 || pr.first==161 || pr.first==163 || pr.first==165 || pr.first==167)
      return 5; // long inter-ganged pixel
    else if (pr.first==153 || pr.first==155 || pr.first==157 || pr.first==159 || pr.first==160 || pr.first==162 || pr.first==164 || pr.first==166)
      return 3; // long ganged pixel
    else
      return 1; // long 
  }
}

bool IsNormalType(int row,int col){  
  return PixelType(row,col) == 0 ? true : false;
}

void ParseCalLine(std::string& line, int& row, int& col,int& tot,int& vcal,int& dummy) {
  
  std::istringstream iss(line);
  //splitted by whitespace
  int dummy1,dummy2;
  iss >> dummy1;
  iss >> row;
  iss >> col;
  iss >> vcal;
  iss >> dummy2;
  iss >> tot;
  iss >> dummy;

  if(debug) {
    cout << "Line "<<line<<" splitted to: " 
	 << dummy1 << " , "
	 << row << " , "
	 << col << " , "
	 << vcal << " , "
	 << dummy2 << " , "
	 << tot << " , "
	 << dummy << endl;
  }
  
} //ParseCalLine

void processFile(const std::string& filename, double cin, std::string name) {
  
  std::map<int,TH1F*> vcal_htot;   
  pixels pixs;
  std::ifstream ifs(filename.c_str());  
  TH2F* h_usedpixels = new TH2F(TString::Format("h_usedpixels_%s",name.c_str()),TString::Format("h_usedpixels_%s",name.c_str()),18,0,18,160,0,160);

  int lines=0;
  while(ifs.good()) {
    std::string line;
    std::getline(ifs,line);
    //int dummy1;
    int row;
    int col;
    int vcal;
    //int dummy2;
    int tot;
    int count;
    ParseCalLine(line, row, col, tot, vcal, count);    
    std::map<int,TH1F*>::iterator ifound = vcal_htot.find(vcal);
    if(debug) cout<<"row="<<row<<" col="<<col<<endl;

    
    if(lines%10000==0) cout<<"Processing line " << lines++ <<endl;
    lines++;
    
    pixel* pix = pixs.FindPixel(row,col);
    if(pix==NULL) {
      pix = pixs.Add(row,col);
      if(debug) cout<<"Added pixel ["<<row<<","<<col<<"]"<<endl;
    }
    pix->Fill(vcal,tot,count);
    
    if(IsMasked(row,col)) {
      if(debug) cout<<"masked!"<<endl;      
      continue;
    }
    
    if(!IsNormalType(row,col)) {
      if(debug) cout<<"is STRANGE ("<<PixelType(row,col)<<")"<<endl;
      continue;
    }    
    if(debug) cout<<"is normal ("<<PixelType(row,col)<<")"<<endl;
    
    if(vcal==830) {
      for(int i=0;i!=count;i++) h_usedpixels->Fill(col,row);      
    }
    
    if(ifound != vcal_htot.end()) {
      if(debug) cout<<"found histo at "<<ifound->second<<endl;
      for(int i=0;i!=count;i++) (ifound->second)->Fill(tot);      
    }else {
      vcal_htot[vcal] = new TH1F(TString::Format("h_tot_vcal%i_%s",vcal,name.c_str()),TString::Format("ToT all pixels vcal=%i %s;ToT;Entries",vcal,name.c_str()),170,-10,160);
      if(debug) cout<<"created histo at "<<vcal_htot[vcal]<<endl;
      vcal_htot[vcal]->Fill(tot);
    }
    
  }//ifs.good

  pixs.MakeCalAll();

  if(debug)
    pixs.PlotCalAll();
  
  pixs.SaveAll(outfile);

  TDirectory* dir = gDirectory;
  outfile->cd();
  for(std::map<int,TH1F*>::iterator i = vcal_htot.begin();i!=vcal_htot.end();++i) {
    i->second->Write();
  } //i
  dir->cd();
  
  //plot a few pixels as example
  for(int irow =25;irow!=31;++irow) {
    pixel* pix = pixs.FindPixel(irow,6);
    if(pix!=NULL) {
      pix->PlotTot(110);
      pix->PlotTot(270);
      pix->PlotTot(430);
      pix->PlotTot(630);
      pix->PlotTot(950);
      pix->PlotCal();
    }
  }
  
  //plot the distribution of the parameters
  TH1F* hp0 = new TH1F("hp0","hp0;parameter p_{0}",100,500,1500);
  TH1F* hp1 = new TH1F("hp1","hp1;parameter p_{1}",100,-1900,-2100);
  TH1F* hp2 = new TH1F("hp2","hp2;parameter p_{2}",100,150000,320000);
  for(std::vector<pixel>::iterator i=pixs.v.begin();i!=pixs.v.end();++i) {
    TF1* f = i->tot_cal_fit;
    if(f!=NULL) {
      hp0->Fill(f->GetParameter(0));
      hp1->Fill(f->GetParameter(1));
      hp2->Fill(f->GetParameter(2));
    }
  }

  TCanvas* c_ps = new TCanvas("c_ps","c_ps",10,20,700,500);
  c_ps->Divide(2,2);
  c_ps->cd(1);
  hp0->Draw();
  c_ps->cd(2);
  hp1->Draw();
  c_ps->cd(3);
  hp2->Draw();
  c_ps->SaveAs(TString::Format("totcal_params_%s.png",name.c_str()));
  delete hp0,hp1,hp2;
  delete c_ps;
  
  
  
  TCanvas* c_pixels = new TCanvas(TString::Format("c_pixels_%s",name.c_str()),TString::Format("c_pixels_%s",name.c_str()),10,20,700,500);
  h_usedpixels->Draw("COLZ");  
  c_pixels->SaveAs(TString::Format("pixels_%s.png",name.c_str()));
  delete h_usedpixels;
  delete c_pixels;
  
  TCanvas* c_all = new TCanvas(TString::Format("c_all_tot_vcal_%s",name.c_str()),TString::Format("c_all_tot_vcal_%s",name.c_str()),10,20,700,500);
  TGraph* gr_entries = new TGraph((int)vcal_htot.size());
  TGraphErrors* gr_tot = new TGraphErrors((int)vcal_htot.size());
  std::map<int,TH1F*>::iterator i  = vcal_htot.begin(); 
  std::map<int,TH1F*>::iterator iE = vcal_htot.end(); 
  int idummy=0;
  for(;i!=iE;++i) {    
    TH1F* h = i->second; 
    if(name=="168"||name=="172") h->Rebin();
    
    bool fit_ok = false;
    while(!fit_ok) {
      h->Fit("gaus","Q");
      TF1* f_tmp = h->GetFunction("gaus");
      //it should be a nice gaussian shape -> other wise rebin
      if(fabs(h->GetMean()-f_tmp->GetParameter("Mean")) < 2.0) {
	fit_ok=true;
      } else {
	//do something
	//rebin
	h->Rebin();
      }
    }
    

    double Q = cin*TMath::Power(10,-15)*i->first*TMath::Power(10,-3);
    double Qe = Q/(1.6*TMath::Power(10,-19));
    //Qe/=1000.0;
    //if(i->first > 300) continue;
    TCanvas* c = new TCanvas(TString::Format("c_tot_vcal%i_%s",i->first,name.c_str()),TString::Format("c_tot_vcal%i_%s",i->first,name.c_str()),10,20,700,500);
    h->Draw();
    
    myText(0.7,0.88,kBlack,TString::Format("VCal=%imV",i->first));
    myText(0.7,0.81,kBlack,TString::Format("Q_{inj}=%.0fe^{-}",Qe));
    myText(0.15,0.88,kBlack,TString::Format("C_{in}=%.1ffF",cin));
    myText(0.75,0.22,kBlack,TString::Format("%s",name.c_str()));
    
    TF1* fg = h->GetFunction("gaus");
    if(fg==0) {cout<<"no function found!"<<endl; exit(0);}
    double C = fg->GetParameter(0);
    double m = fg->GetParameter(1);
    double m_error = fg->GetParError(1);
    double s = fg->GetParameter(2);
    //double s_error = fg->GetParameter(1);
    myText(0.15,0.83,kBlack,TString::Format("C=%.0f,m=%.1f,#sigma=%.1f",C,m,s));
    c->SaveAs(TString::Format("tot_vcal%i_%s.png",i->first,name.c_str()));
    //c->SaveAs(TString::Format("tot_vcal_%i.png",i->first));
    delete c;

    c_all->cd();
    h->SetLineColor(idummy+1);
    if(idummy==0) { h->Draw(); }
    else { h->Draw("same"); }
    idummy++;
    
    gr_entries->SetPoint(idummy-1,i->first,h->GetEntries());
    
    gr_tot->SetPoint(idummy-1,Qe,m);
    gr_tot->SetPointError(idummy-1,Qe*0.01,m_error);
    
    
  }

  c_all->SaveAs(TString::Format("tot_vcal_all_%s.png",name.c_str()));
  delete c_all;

  TCanvas* c_entries = new TCanvas(TString::Format("c_entries_tot_vcal_%s",name.c_str()),TString::Format("c_entries_tot_vcal_%s",name.c_str()),10,20,700,500);
  gr_entries->SetTitle("Injections for each VCal;VCal [mV];Injections");
  gr_entries->Draw("ALP");
  myText(0.75,0.22,kBlack,TString::Format("%s",name.c_str()));
  c_entries->SaveAs(TString::Format("entries_vcal_%s.png",name.c_str()));
  
  TCanvas* c_tot = new TCanvas(TString::Format("c_tot_vcal_%s",name.c_str()),TString::Format("c_tot_vcal_%s",name.c_str()),10,20,700,500);
  gr_tot->SetTitle(";Injected charge [e^{-}];ToT");
  gr_tot->SetMarkerSize(0.9);

  //fit with function
  //use straight line now but get a better function later(?)
  bool polfit = false;
  TF1* f = 0;
  if(polfit) {
    gr_tot->Fit("pol2");
    gr_tot->Draw("AP");
    f = gr_tot->GetFunction("pol2");
    myText(0.75,0.22,kBlack,TString::Format("%s",name.c_str()));
    myBoxMultiLineText(0.15,0.83,TString::Format("a=%.1f|b=%.1f|c=%.1f",
						 f->GetParameter(0),
						 f->GetParameter(1),
						 f->GetParameter(2)));
  }
  else {
    f = new TF1("ftot",totqfnc,4000,45000,3);
    f->SetParameter(0,1000);
    f->SetParameter(1,-1000);
    f->SetParameter(2,500000);    
    gr_tot->Fit("ftot","QR");
    gr_tot->Draw("AP");
    myText(0.75,0.22,kBlack,TString::Format("%s",name.c_str()));
    myBoxMultiLineText(0.15,0.83,TString::Format("p_{0}=%.1f|p_{1}=%.1f|p_{2}=%.1f",
						 f->GetParameter(0),
						 f->GetParameter(1),
						 f->GetParameter(2)));
    
  }
  c_tot->SaveAs(TString::Format("q_vs_tot_%s.png",name.c_str()));
  if(debug) cout << "gpad fillcolor="<<gPad->GetFillColor()<< " and tcanvas fill color="<<   c_tot->GetFillColor() <<endl;

  outfile->cd();
  f->Write(TString::Format("ftotToQ_%s",name.c_str()));
  gr_tot->Write(TString::Format("gr_tot_%s",name.c_str()));
  
  
  //delete f;
  delete c_entries;
  delete gr_entries;
  delete gr_tot;
  delete c_tot;
  
  for(std::map<int,TH1F*>::iterator i=vcal_htot.begin();i!=vcal_htot.end();++i) {
    delete i->second;
  }
  
  return;
  
} //processFile




int main(int argc, char* argv[]) {

  std::string filename,maskfile;
  Cinj=-1;
  if(argc != 4 && argc !=5 && argc !=6) {
    cout << argv[0] 
	 << " requires 3, 4 or 5 arguments: \n"
	 << "\tEx: "<<argv[0]<<" raw_tot_cal_file.tot Cinj sensorname [mask file] [debug]" << endl;
    exit(-1);      
  } else {
    filename = argv[1];
    Cinj = atof(argv[2]);
    name = argv[3];
    if(argc==5) {
      std::string s = argv[4];
      if(s.find("debug")==std::string::npos) maskfile=argv[4];
      else debug = 1;
    } else {
      maskfile=argv[4];
      debug=1;
    }
  }

  TDirectory* dir = gDirectory;
  outfile = new TFile(TString::Format("totToQ_%s.root",name.c_str()),"RECREATE");
  dir->cd();
  
  cout << "Executing:" 
       << "\nName\t"<<name
       << "\nToT cal file\t"<<filename
       << "\nCinj\t"<<Cinj<<" fF"
       << "\nMask file\t"<<maskfile
       << "\nDebug\t"<<debug
       <<"\n\nOutput is saved to: "<< outfile->GetName()<<endl;
  
  rootlogon();
  
  //add masks
  if(maskfile!="") addMask(maskfile);
  
  //   //open the files to be used
  //   for(int isensor=0;isensor!=1;++isensor) {
  //     std::string filename = path+fileNames[isensor];
  processFile(filename,Cinj,name);     
  //  } //isensor
  
  outfile->Close();  
  return 1;
}






