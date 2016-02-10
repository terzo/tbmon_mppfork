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

#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"

/** \brief description of Botho analysis
*/

class Botho: public TbAnalysis {

private:

	////////////////////////////
	// DATA ////////////////////
	////////////////////////////
	
	double pitchX;	///< pitchX of the DUT
	double pitchY;	///< pitchY of the DUT
	double epitchX;	///< pitchX of edge pixels of DUT
	double epitchY;	///< pitchY of edge pixels of DUT
	
	double n_cols;	///< number of columns of DUT
	double n_rows;	///< number of rows of the DUT
	
	int skip_cols;	///< number of columns to be skipped in analysis
	int skip_rows;	///< number of rows to be skipped in analysis
	
	int lv1Min;	///< lvl1 minimum limit of DUT
	int lv1Max;	///< lvl1 maximum limit of DUT
	
	bool isFEI4;	///< true if FE-I4, false otherwise, needs to be set
	bool isFEI3;	///< true if FE-I3, false otherwise, needs to be set
	
	event::eventflag_t (Botho::*track_check_central_region)(const Event &);	///< pointer to cut function
	
	
	////////////////////////////
	// CONSTANTS ///////////////
	////////////////////////////
	
	
/*	// region of interest (ROI) in pixels
	// 1cmx1cm around guessed beam center of (53,195)
	// for 250umx50um pixels this is 40x200 pixels
	static const int ROI_x1 = 33;
	static const int ROI_y1 = 95;
	static const int ROI_x2 = 72;
	static const int ROI_y2 = 294;
*/	
	// define square coordinate array
	int ROI_square[4];	// = {ROI_x1, ROI_y1, ROI_x2, ROI_y2};

	
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
	TH2D* h_hitmap_matched_pos;	///< positions of matched hits
	TH1D* h_profileX;	///< profile of hits along columns
	TH1D* h_profileY;	///< profile of hits along rows
	TH1D* h_profileX_pixel;	///< profile of hits along columns folded into pixel
	TH1D* h_profileY_pixel;	///< profile of hits along rows folded into pixel
	// beam profiles with 10, 5, or 2 bins per col/row pitch
	TH2D* h_beamProfileFEI3;	///< beam profile
	TH2D* h_beamProfileFEI4;
//	TH2D* h_beamProfile2;
	// angle profiles
	TH2D* h_angleProfile;	///< shows a map of average angle per position
	TH2D* h_angleProfile_numEvents;	///< auxiliary histogram for the angle profi
	// residual histograms
	TH1D *h_resX;
	TH1D *h_resY;
	// pixel profiles with cuts and residuals
	TH1D *h_profileX_pixel_rescuts;
	TH1D *h_profileY_pixel_rescuts;
	TH1D *h_profileX_pixel_goodtrackcuts;
	TH1D *h_profileY_pixel_goodtrackcuts;
	TH1D *h_resX_CS1_unweighted;
	TH1D *h_resY_CS1_unweighted;
	TH1D *h_resX_CS2_unweighted;
	TH1D *h_resY_CS2_unweighted;
	TH1D *h_resX_CSall_MaxTotCell;
	TH1D *h_resY_CSall_MaxTotCell;
	TH1D *h_resX_CSall_CXonly_MaxTotCell;
	TH1D *h_resY_CSall_CYonly_MaxTotCell;
	TH1D *h_resX_CS2x2max_MaxTotCell;
	TH1D *h_resY_CS2x2max_MaxTotCell;
	// module profiles with cuts
	event::eventflag_t (Botho::*track_geometry_and_mask)(const Event &);	///< pointer to cut function
	TH1D *h_profileX_module_goodtrackcuts;
	TH1D *h_profileY_module_goodtrackcuts;
	TH1D *h_profileX_module_rescuts;
	TH1D *h_profileY_module_rescuts;
	// pixel angle profiles
	TH1D *h_angleProfileX_rescuts_avg;
	TH1D *h_angleProfileX_rescuts_hits;
	TH1D *h_angleProfileY_rescuts_avg;
	TH1D *h_angleProfileY_rescuts_hits;
	
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
	TH1D* h_clusize_x;
	TH1D* h_clusize_x_CSx;
	TH1D* h_clusize_y;
	TH1D* h_clusize_y_CSy;
	TH1D* h_cluster_multiplicity;	///< number of clusters in one event for Good events (event.fTrack == event::kGood)
	TH1D* h_cluster_multiplicity_matched;	///< number of clusters in matched events
	// cluster track maps
	static const int num_of_cs_track_histos = 6;
	TH2D* h_track_cs[num_of_cs_track_histos];	///< track of cluster size 1-5 and more
//	TH2D* h_track_1x1um_cs[num_of_cs_track_histos];
	// cluster track maps projected into one pixel; use CUT3: event.fRegion == event::kGood
	TH2D* h_pixel_track_cs_all;	///< pixelmap of all hits
	TH2D* h_pixel_track_1x1um_cs_all;	///< pixelmap 1x1um bins of all hits
	TH2D* h_pixel_track_cs[num_of_cs_track_histos];	///< track of cluster size 1-5 and more
	TH2D* h_pixel_track_1x1um_cs[num_of_cs_track_histos];	//< track of cluster size 1-5 and more with 1x1um bins
	
	// map resolution (for rebinning) in um
	static const int cs_pix_maps_res_x = 5;
	static const int cs_pix_maps_res_y = 5;
	TH2D* h_pixel_track_cs_all_twopix;	///< pixelmap of all hits over two pixel distance
	TH2D* h_pixel_track_cs_twopix[num_of_cs_track_histos];
	
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
	// ToT histograms
	TH1D *h_tot_all;
	TH1D *h_tot_cs1;
	TH1D *h_tot_cs2;
	// ToT module map
	TH2D *h_module_tot_hits;
	TH2D *h_module_tot_avg;
	// ToT module map cluster size 1 hits only
	TH2D *h_module_tot_hits_cs1;
	TH2D *h_module_tot_avg_cs1;
	// ToT module map cluster size 2 hits only
	TH2D *h_module_tot_hits_cs2;
	TH2D *h_module_tot_avg_cs2;
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
	TH1D *h_tot_all_selected[num_of_tot_selected_divisions+2];
	TH1D *h_tot_cs1_selected[num_of_tot_selected_divisions+2];
	TH1D *h_tot_cs2_selected[num_of_tot_selected_divisions+2];
	// square selection
	TH2D *h_module_totmap_square_selected [2];
	TH2D *h_pixel_tot_hits_square_selected[2];
	TH2D *h_pixel_tot_avg_square_selected [2];
	TH1D *h_tot_all_square_selected[2];
	TH1D *h_tot_cs1_square_selected[2];
	TH1D *h_tot_cs2_square_selected[2];
	// ToT pixel map
	TH2D *h_pixel_tot_avg;
	TH2D *h_pixel_tot_hits;
	// ToT pixel maps with cluster size cuts
	TH2D *h_pixel_tot_hits_CS[4];
	TH2D *h_pixel_tot_avg_CS[4];
	// alternative ToT pixel map
	TH2D *h_pixel_tot_avg_alt;
	TH2D *h_pixel_tot_hits_alt;
	// ToT 3D scatter plot with 5x5um binning
	TH3D *h_pixel_tot_3d;
	// high cluster size map of matched clusters
	TH2D *h_cs09_map;
	TH2D *h_cs10_map;
	TH2D *h_cs11_map;
	/// ---- COMMAND LINE SWITCHES ---- ///
	bool create_totmap_file;	// determines if totmap file is created
	bool do_totmap_selection;	// does analysis depending no totmap file values
	bool do_square_selection;	// does analysis with square area selection
	
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
	TH2D *h_trackPixelMap_square_selected;
	TH2D *h_hitPixelMap_square_selected;
	// efficiency2.cc rebuild data and histograms
	int eff2cc_ntracks;
	int eff2cc_nhits;
	int eff2cc_anyhit;
	TH2D *h_eff2cc_trackMap;
	TH2D *h_eff2cc_hitMap;
	TH2D *h_eff2cc_anyhitMap;
	// efficiency2.cc pixel maps
	TH2D *h_eff2cc_trackPixelMap;
	TH2D *h_eff2cc_hitPixelMap;
	TH2D *h_eff2cc_trackPixelMap_square_selected;
	TH2D *h_eff2cc_hitPixelMap_square_selected;
	/// ---- COMMAND LINE SWITCHES ---- ///
	bool do_eff_square_selection;	// does analysis with square area selection


	///////////////////////////////////////
	/// EDGE EFFICIENCY ANALYSIS ----------
	///////////////////////////////////////
	// initialize histograms etc
	void init_edgeefficiency_analysis(TbConfig &config);
	// do analysis of event
	void edgeefficiency_analysis(const TbConfig &config, const Event &event);
	// save histograms and free memory
	void finalize_edgeefficiency_analysis(const TbConfig &config);
	/// ---- DATA ---- ///
	int ee_ntracks;
	int ee_nTracksEdge;
	int ee_nHitsEdge;
	TH1D *h_hitOutLeft;
	TH1D *h_trackOutLeft;
	TH1D *h_effOutLeft;
	TH1D *h_hitOutRight;
	TH1D *h_trackOutRight;
	TH1D *h_effOutRight;
	
	TH2D* h_EEAnyHitSensorMap;
	TH2D* h_EEHitSensorMap;
	TH2D* h_EETrackSensorMap;
	TH2D* h_EEHitPixelMap;
	TH2D* h_EEHitLeftPixelMap;
	TH2D* h_EEHitRightPixelMap;
	TH2D* h_EETrackPixelMap;
	TH2D* h_EETrackLeftPixelMap;
	TH2D* h_EETrackRightPixelMap;
	
	TH1D *h_bothoEE_hitOutLeft;
	TH1D *h_bothoEE_trackOutLeft;
	TH1D *h_bothoEE_effOutLeft;
	TH1D *h_bothoEE_hitOutRight;
	TH1D *h_bothoEE_trackOutRight;
	TH1D *h_bothoEE_effOutRight;
	/// ---- FUNCTIONS ---- ///
	void calc_efficiency(double, double, double &, double &);
	void drawToFileEdgeEff(const TbConfig& config,
		const char* analysisName, const char* histoName, TObject* histo,
		double fitLeft, double fitRight) const;
	
private:
	// private management functions
	
	// own saveToFile function creating subdirectories
	void saveToFile(const TbConfig &config, const char *analysisName,
			const char *subfolderName, const char *histoName, TNamed *histo);
	// function checking square
	bool in_square(int *square, int x, int y);
	// function drawing pixel with correct aspect ratio
	void drawPixelNicely(const TbConfig&, const char*, TH1*);
	
	
	// geometric cuts for different modules (finer than fTrackCentraRegion)
	event::eventflag_t check_fTrackCentralRegion(const Event &event) {
		return event.fTrackCentralRegion; }
	event::eventflag_t check_fTrackRegion(const Event &event) {
		return event.fTrackRegion; }
	event::eventflag_t check_central_region_VTT_FEI3_TBNov2013(const Event &event);
	event::eventflag_t check_central_region_VTT_FEI3_TBNov2013_2(const Event &event);
	event::eventflag_t check_fTrackRegion_geoid100_referenceDUT(const Event &event);
	event::eventflag_t check_fTrackRegion_geoid198to204_DUT(const Event &event);

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
	
};

#endif //BOTHO_H
