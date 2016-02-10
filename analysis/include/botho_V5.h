/**  \file botho.h
	\brief header file for Botho's analysis
	
	Botho's analysis
	written by Botho Paschen bpaschen@mppmu.mpg.de
*/
#ifndef BOTHO_H
#define BOTHO_H

//TBmon headers
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

//standard headers
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

//root headers
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TPaletteAxis.h"

/** \brief description of Botho analysis
*/

class Botho: public TbAnalysis {

private:

	////////////////////////////
	// DATA ///////////////////
	////////////////////////////
	
	double pitchX;	///< pitchX of the DUT
	double pitchY;	///< pitchY of the DUT
	
	double n_cols;	///< number of columns of DUT
	double n_rows;	///< number of rows of the DUT
	
	int lv1Min;	///< lvl1 minimum limit of DUT
	int lv1Max;	///< lvl1 maximum limit of DUT
	
	////////////////////////////
	/// HISTOGRAMS /////////////
	////////////////////////////
	
	TH2D* h_masks;		///< histogram that shows the masks that were read in
					// by DUT->addMasks() into the
					// DUT::masks<std::pair<int, int> > vector
	TH2D* h_lv1excluded;	///< contains all raw hits that get excluded by lv1
						// but are not excluded by the global mask
	
	////////////////////////////
	/// SUBANALYSES ////////////
	////////////////////////////
	
	
	///////////////////////////////////////
	/// HIT ANALYSIS ----------------------
	///////////////////////////////////////
	// initialize histograms etc
	void init_hit_analysis(TbConfig &config);
	// do analysis of event
	void hit_analysis(const TbConfig &config, const Event &event);
	// save histograms and free memory
	void finalize_hit_analysis(const TbConfig &config);
	/// ---- DATA ---- ///
	TH1D* h_hit_multiplicity;	///< number of hits per event
	TH1D* h_hit_multiplicity_matched;	///< number of hits in matched events
//	this was just for checking, whether all events in the analysis for one DUT
//	are in the correct DUT. iden is the DUT ID number
	TH1D* h_hit_iden;	///< hit iden
	TH2D* h_hitmap_all;	///< column and row position of all hits
	TH2D* h_hitmap_all_check;	///< check of hitmap loops
	TH2D* h_hitmap_all_pos;	///< more accurate position of all hits
	TH2D* h_hitmap_raw_pos;	///< raw hitmap
	TH2D* h_rawHits;	///< all raw hits of events
	TH1D* h_profileX;	///< profile of hits along columns
	TH1D* h_profileY;	///< profile of hits along rows
	TH1D* h_profileX_pixel;	///< profile of hits along columns folded into pixel
	TH1D* h_profileY_pixel;	///< profile of hits along rows folded into pixel
	// beam profiles with 10, 5, or 2 bins per col/row pitch
	TH2D* h_beamProfile10;	///< beam profile
	TH2D* h_beamProfile5;
	TH2D* h_beamProfile2;
	// angle profiles
	TH2D* h_angleProfile;	///< shows a map of average angle per position
	TH2D* h_angleProfile_numEvents;	///< auxiliary histogram for the angle profi
	// residual histograms
	TH1D *h_resX;
	TH1D *h_resY;
	
	
	///////////////////////////////////////
	/// CLUSTER ANLAYSIS ------------------
	///////////////////////////////////////
	// initialize histograms etc
	void init_cluster_analysis(TbConfig &config);
	// do analysis of event
	void cluster_analysis(const TbConfig &config, const Event &event);
	// save histograms and free memory
	void finalize_cluster_analysis(const TbConfig &config);
	/// ---- DATA ---- ///
	// cluster histograms -----------------------------------------------------
	// all clusters of events with flag event.fTrack == event::kGood
	TH1D* h_clusize_all;
	TH1D* h_clusize_matched;	///< matched clusters with flag event.fTrack == event::kGood
	TH1D* h_clusize_unmatched;	///< all unmatched clusters with flag event.fTrack == event::kGood
	TH1D* h_clusize_matched_good_region;	///< Total cluster size distribution of all clusters which passed the following cuts: central region and matched (taken from Residuals Analysis)
	TH1D* h_cluster_multiplicity;	///< number of clusters in one event for Good events (event.fTrack == event::kGood)
	TH1D* h_cluster_multiplicity_matched;	///< number of clusters in matched events
	// cluster track maps
	static const int num_of_cs_track_histos = 6;
	TH2D* h_track_cs[num_of_cs_track_histos];	///< track of cluster size 1-5 and more
	// cluster track maps projected into one pixel; use CUT3: event.fRegion == event::kGood
	TH2D* h_pixel_track_cs[num_of_cs_track_histos];	///< track of cluster size 1-5 and more
	
	
	///////////////////////////////////////
	/// TOT - CHARGE ANALYSIS -------------
	///////////////////////////////////////
	// initialize histograms etc
	void init_tot_analysis(TbConfig &config);
	// do analysis of event
	void tot_analysis(const TbConfig &config, const Event &event);
	// save histograms and free memory
	void finalize_tot_analysis(const TbConfig &config);
	/// ---- DATA ---- ///
	// ToT module map
	TH2D *h_module_tot_hits;
	TH2D *h_module_tot_avg;
	// optional tot selection map
	static const double tot_selected_lowerbound=1;
	static const double tot_selected_upperbound=7;
	static const int num_of_tot_selected_divisions = 1;
	// set span of tot for selection
	double tot_sel_span;	//= std::fabs(tot_selected_upperbound - tot_selected_lowerbound);
	// calculate size of one "bin"
	double tot_sel_division;	//= tot_sel_span/num_of_tot_selected_divisions;
	TH2D *h_module_tot_selection;
	// selected maps include 1 underflow and 1 overflow map
	TH2D *h_module_totmap_selected [num_of_tot_selected_divisions+2];
	TH2D *h_pixel_tot_hits_selected[num_of_tot_selected_divisions+2];
	TH2D *h_pixel_tot_avg_selected [num_of_tot_selected_divisions+2];
	// ToT pixel map
	TH2D *h_pixel_tot_avg;
	TH2D *h_pixel_tot_hits;
	// alternative ToT pixel map
	TH2D *h_pixel_tot_avg_alt;
	TH2D *h_pixel_tot_hits_alt;
	// ToT 3D scatter plot with 5x5um binning
	TH3D *h_pixel_tot_3d;
	/// ---- COMMAND LINE SWITCHES ---- ///
	bool create_totmap_file;	// determines if totmap file is created
	bool do_totmap_selection;	// does analysis depending no totmap file values
	
	
	///////////////////////////////////////
	/// CUT ANALYSIS -----------------------
	///////////////////////////////////////
	// initialize histograms etc
	void init_cut_analysis(TbConfig &config);
	// do analysis of event
	void cut_analysis(const TbConfig &config, const Event &event);
	// save histograms and free memory
	void finalize_cut_analysis(const TbConfig &config);
	/// ---- DATA ---- ///
	TH1D *h_cuts_single;
	TH1D *h_cuts_combined;
	TH1D *h_cuts_botho_single;
	TH1D *h_cuts_botho_combined;
	
	
	///////////////////////////////////////
	/// EFFICIENCY ANALYSIS ----------------
	///////////////////////////////////////
	// initialize histograms etc
	void init_efficiency_analysis(TbConfig &config);
	// do analysis of event
	void efficiency_analysis(const TbConfig &config, const Event &event);
	// save histograms and free memory
	void finalize_efficiency_analysis(const TbConfig &config);
	/// ---- DATA ---- ///
	int n_tracks_raw;	///< number of raw tracks
	int n_tracks_good;	///< number of good tracks (pass CUT 1)
	int n_hits;	///< number of matched hits (for efficiency)
	// efficiency module histograms
	TH2D *h_trackModuleMap;
	TH2D *h_hitModuleMap;
	// efficiency pixel histograms (copied from efficiency2 analysis)
	TH2D *h_trackPixelMap;
	TH2D *h_hitPixelMap;
	// efficiency.cc rebuild data and histograms
	int effcc_ntracks;
	int effcc_nhits;
	int effcc_anyhit;
	TH2D *h_effcc_trackMap;
	TH2D *h_effcc_hitMap;
	TH2D *h_effcc_anyHitMap;

private:
	// private management functions
	
	// own saveToFile function creating subdirectories
	void saveToFile(const TbConfig &config, const char *analysisName,
			const char *subfolderName, const char *histoName, TNamed *histo);
	// function drawing pixel with correct aspect ratio
	void drawPixelNicely(const TbConfig&, const char*, TH1*);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
	
};

#endif //BOTHO_H
