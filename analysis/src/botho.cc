/**  \file botho.cc
	\brief source file for Botho's analysis
	
	Botho's analysis
	written by Botho Paschen bpaschen@mppmu.mpg.de
*/
#include "botho.h"	// include obvious header for the analysis

///////////////////////////////////////////////
/// INITIALIZATION OF ANALYSIS 
///////////////////////////////////////////////
void Botho::init(TbConfig &config) {
        
	//TH1::SetDefaultSumw2();
	// copy DUT information from the config object to members of this object
	const DUT* dut = config.getDut(this->iden);
	
	// get pitches
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();
	
	// get column and row numbers
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();
	
	// get lvl1 limits
	lv1Min = dut->lv1Min;
	lv1Max = dut->lv1Max;
	
	// set isFEI3 and isFEI4 values
	// (this has to be copied to every analysis using this)
	// check for FE-I4
	if( pitchX==250 && n_cols == 80 ) {
		isFEI4 = true;
		TBALOG(kINFO) << "module recognized as FE-I4" << endl;
	}
	else isFEI4 = false;
	// check for FE-I3
	if( pitchX==400 && n_cols == 18 ) {
		isFEI3 = true;
		TBALOG(kINFO) << "module recognized as FE-I3" << endl;
	}
	else isFEI3 = false;
	if(isFEI3 && isFEI4)
		TBALOG(kERROR) << "module is FE-I3 AND FE-I4???" << endl;

	// define central region cut by setting function pointer
	// (this has to be copied to every analysis using this)
	if(strcmp(config.tbslot, "eudetPPSNov2013FEI34") == 0
		|| strcmp(config.tbslot, "eudetPPSSep2012FEI3") == 0) {
		if(isFEI3) {
			track_check_central_region = &Botho::check_central_region_VTT_FEI3_TBNov2013_2;  // changed 2014-10-28
			TBALOG(kINFO) << "initiated special geometric FE-I3 cut"
					<< " for November 2013 / September 2012 TB modules" << endl;
		}
		else
			track_check_central_region = &Botho::check_fTrackCentralRegion;
	}
	else
		track_check_central_region = &Botho::check_fTrackCentralRegion;
	
	// SPECIAL region of interest (ROI) stuff for highly irradiated CIS
	// region of interest (ROI) in pixels
	// 1cmx1cm around guessed beam center of (53,195)
	// for 250umx50um pixels this is 40x200 pixels
	const int ROI_x1 = 33;
	const int ROI_y1 = 95;
	const int ROI_x2 = 72;
	const int ROI_y2 = 294;
		
	// define square coordinate array
	// ROI_square[4] = {ROI_x1, ROI_y1, ROI_x2, ROI_y2};
	ROI_square[0] = ROI_x1;
	ROI_square[1] = ROI_y1;
	ROI_square[2] = ROI_x2;
	ROI_square[3] = ROI_y2;

	
	///////////////////////////
	// INITIALIZE HISTOGRAMS -------------------------------------------------
	/// define and initialize histograms
	
	/// lv1 excluded map
	h_lv1excluded = new TH2D("", ";col;row", n_cols, 0, n_cols, n_rows, 0, n_rows);
	
	/// mask map		
	h_masks = new TH2D("", ";col;row", n_cols, 0, n_cols, n_rows, 0, n_rows);
	// routine copied from maskmaker.cc and customized
	for(int col = 0; col < n_cols; col++) {
		for(int row = 0; row < n_rows; row++) {
			for (vector<pair<int, int> >::const_iterator mm = dut->masks.begin();
						mm != dut->masks.end(); mm++) {
				// corrected col and row
				if ((*mm).first == col and (*mm).second == row) {
					h_masks->Fill(col, row);
					break;
				}
			}
		}
	}
	
	/// initialize hit analysis
	init_hit_analysis(config);
			
	/// initialize cluster analysis
	init_cluster_analysis(config);
	
	/// initialize tot analysis
	init_tot_analysis(config);
	
	/// initialize efficiency analysis
	init_efficiency_analysis(config);
			
	/// initialize cut analysis
	init_cut_analysis(config);
	
	/// initialize edge efficiency analysis
	init_edgeefficiency_analysis(config);
}


////////////////////////////////////////////////////////////////////////////////
/// EVENT PROCESSING ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Botho::event(const TbConfig &config, const Event &event) {
	/// do hit analysis
	hit_analysis(config, event);

	/// do cluster analysis
	cluster_analysis(config, event);

	/// do tot analysis
	tot_analysis(config, event);
	
	/// do cut analysis
	cut_analysis(config, event);

	/// do efficiency analysis
	efficiency_analysis(config, event);
	
	/// do edge efficiency analysis
	edgeefficiency_analysis(config, event);
	

	////////////////////////
	//// RAW DATA ANALYSIS ----------------------------------------------------
	////////////////////////
	
	// raw hits
	for(vector<PllHit *>::const_iterator it=event.rawHits.begin();
			it != event.rawHits.end(); it++) {
		int col = (*it)->col;
		int row = (*it)->row;
		int lv1 = (*it)->lv1;
		// fill h_lv1excluded if raw hit is not excluded by mask, but by lv1
		if((h_masks->GetCellContent(col, row) == 0)
				&& ((lv1 < lv1Min) || (lv1 > lv1Max))) {
			h_lv1excluded->Fill(col, row);
		}
	}
	

	////////////////////////////////////////////
	/// CUT 1:  EVENT.FTRACK == EVENT::KGOOD ///
	////////////////////////////////////////////
	
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {  // if not event.fTrack == event::kGood
		return;                         // leave event
	}
	TBALOG(kDEBUG3) << "Event passed cut 1: \"event.fTrack == event::kGood\"."
	                << endl;      // debug output message
	
	/////////////////////
	//// CUT 1 ANALYSIS -------------------------------------------------------
	/////////////////////
	   
	
	//////////////////////////////////
	/// CUT 2: GET MATCHED CLUSTER ///
	//////////////////////////////////
	
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
	                       // of matched cluster in vector of clusters of event
	                       
	if (matched == -1) {  // in case there is no matched cluster
		// debug output message
		TBALOG(kDEBUG2) << "No matched clusters --> leaving event." << endl;
		                                             
		return;	// leave event
	}
	// debug output message
	TBALOG(kDEBUG3) << "Event passed cut 2: matched cluster. Additional info:"
	                << "Found " << matched << " matched clusters." << endl;
		
	// get reference to the matched cluster
	const vector<PllHit*>& matched_cluster = event.clusters.at(matched);
	
	
	/////////////////////
	//// CUT 2 ANALYSIS -------------------------------------------------------
	/////////////////////
	
	
	//////////////////////////////////////////////////
	/// CUT 3:  EVENT.FTRACKREGION == EVENT::KGOOD ///
	//////////////////////////////////////////////////
	
	// Track in good region? (central and unmasked) (see event.h)
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	TBALOG(kDEBUG3) << "Event passed cut 3: "
	                << "\"event.fTrackRegion == event::kGood\"." << endl;	// debug output
	
	
	/////////////////////	
	//// CUT 3 ANALYSIS -------------------------------------------------------
	/////////////////////
	
}


////////////////////////////////////////////////////////////////////////////////
/// FINALIZE ANALYSIS AFTER COMPLETE PROCESSING ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Botho::finalize(const TbConfig &config) {

	/////////////////////////////////////////
	/// set viewing ranges for histograms ///
	/////////////////////////////////////////
	h_clusize_all->GetXaxis()->SetRange(1, 6);
	h_clusize_matched->GetXaxis()->SetRange(1, 6);	/// 
	h_clusize_unmatched->GetXaxis()->SetRange(1, 6);
	h_clusize_matched_good_region->GetXaxis()->SetRange(1, 6);
	
	h_hit_multiplicity->GetXaxis()->SetRange(0, 10);
	h_hit_multiplicity_matched->GetXaxis()->SetRange(0, 10);


	/////////////////////////////////////////
	/// print and save histograms ///////////
	/////////////////////////////////////////
	
	// outside lvl1cut
	config.drawAndSave(this->name, "lvl1_excluded", h_lv1excluded);
	
	// mask map
	config.drawAndSave(this->name, "masks", h_masks);
	
	// finalize hit analysis
	finalize_hit_analysis(config);
 	
 	// finalize cluster analysis
 	finalize_cluster_analysis(config);
 	
 	// finalize ToT analysis
 	finalize_tot_analysis(config);
	
	// finalize cuts analysis
	finalize_cut_analysis(config);
	
	// finalize efficiency analysis
	finalize_efficiency_analysis(config);
	
	// finalize edge efficiency analysis
	finalize_edgeefficiency_analysis(config);

	
	////////////////////////////////
	/// free space of histograms ///
	////////////////////////////////
	
	delete h_lv1excluded;
	delete h_masks;
}

////////////////////////////////////////////////////////////////////////////////
/// HIT ANALYSIS            ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// initialize histograms etc for HIT ANALYSIS
void Botho::init_hit_analysis(TbConfig &config) {
	// copy DUT information from the config object to members of this object
	const DUT* dut = config.getDut(this->iden);
	// get pitches
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();	
	// get column and row numbers
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();
	
	// define geometry and mask cut by setting function pointer
	// (this has to be copied to every analysis using this)
	track_geometry_and_mask = &Botho::check_fTrackRegion;
/*	TBALOG( kINFO ) << config.tbslot << endl
				 << this->iden << endl
				 << config.currentRun << endl;
*/
	if(strcmp(config.tbslot, "eudetPPSNov2013FEI34") == 0 && this->iden == 20
			&& strcmp(config.name, "geoid100") == 0) {   // config.currentRun >= 1159 && config.currentRun <= 1187) {
		TBALOG( kINFO ) << "special geometric cut applied for DUT 20 of GEOID100" << std::endl;
		track_geometry_and_mask = &Botho::check_fTrackRegion_geoid100_referenceDUT;
	}
	if(strcmp(config.tbslot, "eudetPPSFeb2014FEI4") == 0 && this->iden == 20
			&& (strcmp(config.name, "geoid198") == 0
			||  strcmp(config.name, "geoid199") == 0
			||  strcmp(config.name, "geoid200") == 0
			||  strcmp(config.name, "geoid201") == 0
			||  strcmp(config.name, "geoid202") == 0
			||  strcmp(config.name, "geoid203") == 0
			||  strcmp(config.name, "geoid204") == 0)) {   // config.currentRun >= 1159 && config.currentRun <= 1187) {
		TBALOG( kINFO ) << "special geometric cut applied for DUT 20 of " << config.name << std::endl;
		track_geometry_and_mask = &Botho::check_fTrackRegion_geoid198to204_DUT;
	}
	
	/// hit multiplicity histograms
	h_hit_multiplicity = new TH1D("", ";number of hits;events",
		100, 0, 100);
	h_hit_multiplicity_matched = new TH1D("", ";number of hits;events",
		100, 0, 100);
		
	/// hit iden histogram
	h_hit_iden = new TH1D("", ";iden;events", 100, 0, 100);
	
	/// hit maps
	int res1 = 10;	// resolution enhancement of the hitmap position with respect to columns and rows
	h_rawHits = new TH2D("", ";column;row",
			n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_profileX = new TH1D("", ";x [#mum];Hits",
			(n_cols+4)*pitchX/10, -2*pitchX, (n_cols+2)*pitchX);
	h_profileY = new TH1D("", ";y [#mum];Hits",
			(n_rows+8)*pitchY/10, -4*pitchY, (n_rows+4)*pitchY);
	h_profileX_pixel = new TH1D("", ";x [#mum];Hits",
			pitchX, 0, pitchX);
	h_profileY_pixel = new TH1D("", ";y [#mum];Hits",
			pitchY, 0, pitchY);
	h_hitmap_all = new TH2D("", ";column;row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
	h_hitmap_all_check = new TH2D("", ";column;row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
	h_hitmap_raw_pos = new TH2D("", ";x [um];y [um]",
			(n_cols+4)*res1, -2*pitchX, (n_cols+2)*pitchX,
			(n_rows+4)*res1, -2*pitchY, (n_rows+2)*pitchY);
	h_hitmap_all_pos = new TH2D("", ";x [um];y [um]",
			(n_cols+4)*res1, -2*pitchX, (n_cols+2)*pitchX,
			(n_rows+4)*res1, -2*pitchY, (n_rows+2)*pitchY);
			
	// matched hitmap position
	h_hitmap_matched_pos = new TH2D("", ";x [#mum]; y [#mum]",
			(n_cols+4)*res1, -2*pitchX, (n_cols+2)*pitchX,
			(n_rows+4)*res1, -2*pitchY, (n_rows+2)*pitchY);
			
	/// beam profiles	
	res1 = 10;	// FE-I3 optimized with 5x5um^2 bins
	h_beamProfileFEI3 = new TH2D("", ";x [um];y [um]",
			(n_cols*3)*res1*8, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1, -1*n_rows*pitchY, (n_rows*2)*pitchY);
//	res1 = 10;	// FE-I4 optimized with 5x5um^2 bins
	h_beamProfileFEI4 = new TH2D("", ";x [um];y [um]",
			(n_cols*3)*res1*5, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1, -1*n_rows*pitchY, (n_rows*2)*pitchY);
/*	res1 = 2;
	h_beamProfile2 = new TH2D("", ";x [um];y [um]",
			(n_cols*3)*res1, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1, -1*n_rows*pitchY, (n_rows*2)*pitchY);
*/
			
	/// angle profile
	res1 = 1;
	h_angleProfile = new TH2D("", ";x [um];y [um];angle",
			(n_cols*3)*res1, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1*0.125, -1*n_rows*pitchY, (n_rows*2)*pitchY);
	h_angleProfile_numEvents = new TH2D("", ";x [um];y [um];angle",
			(n_cols*3)*res1, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1*0.125, -1*n_rows*pitchY, (n_rows*2)*pitchY);
			
	/// residual histograms
	h_resX = new TH1D("", ";x [#mum]", 100, -pitchX * 2, pitchX * 2);
	h_resY = new TH1D("", ";y [#mum]", 100, -pitchY * 2, pitchY * 2);
	
	/// pixel profiles with cuts and resdiuals
	h_profileX_pixel_rescuts = new TH1D("", ";x [#mum]; hits",
			pitchX, 0, pitchX);
	h_profileY_pixel_rescuts = new TH1D("", ";x [#mum]; hits",
			pitchY, 0, pitchY);
	h_profileX_pixel_goodtrackcuts = new TH1D("", ";x [#mum]; hits",
			pitchX, 0, pitchX);
	h_profileY_pixel_goodtrackcuts = new TH1D("", ";x [#mum]; hits",
			pitchY, 0, pitchY);
	h_resX_CS1_unweighted = new TH1D("", ";x [#mum]; hits", 4*pitchX,
			-pitchX * 2, pitchX * 2);
	h_resX_CS2_unweighted = new TH1D("", ";x [#mum]; hits", 4*pitchX,
			-pitchX * 2, pitchX * 2);
	h_resY_CS1_unweighted = new TH1D("", ";y [#mum]; hits", 4*pitchY,
			-pitchY * 2, pitchY * 2);
	h_resY_CS2_unweighted = new TH1D("", ";y [#mum]; hits", 4*pitchY,
			-pitchY * 2, pitchY * 2);
	h_resX_CSall_MaxTotCell = new TH1D("", ";x [#mum];Hits", 4*pitchX,
			-pitchX * 2, pitchX * 2);
	h_resY_CSall_MaxTotCell = new TH1D("", ";y [#mum];Hits", 4*pitchY,
			-pitchY * 2, pitchY * 2);
	h_resX_CSall_CXonly_MaxTotCell = new TH1D("", ";x [#mum];Hits", 4*pitchX,
			-pitchX * 2, pitchX * 2);
	h_resY_CSall_CYonly_MaxTotCell = new TH1D("", ";y [#mum];Hits", 4*pitchY,
			-pitchY * 2, pitchY * 2);
	h_resX_CS2x2max_MaxTotCell = new TH1D("", ";x [#mum];Hits", 4*pitchX,
			-pitchX * 2, pitchX * 2);
	h_resY_CS2x2max_MaxTotCell = new TH1D("", ";y [#mum];Hits", 4*pitchY,
			-pitchY * 2, pitchY * 2);
			
	/// module profiles with cuts
	h_profileX_module_rescuts = new TH1D("", ";x [#mum]; hits",
			n_cols*pitchX/5, 0, n_cols*pitchX);
	h_profileY_module_rescuts = new TH1D("", ";y [#mum]; hits",
			n_rows*pitchY/5, 0, n_rows*pitchY);
	h_profileX_module_goodtrackcuts = new TH1D("", ";x [#mum]; hits",
			n_cols*pitchX/5, 0, n_cols*pitchX);
	h_profileY_module_goodtrackcuts = new TH1D("", ";y [#mum]; hits",
			n_rows*pitchY/5, 0, n_rows*pitchY);
	
	/// angle pixel profiles
	h_angleProfileX_rescuts_avg = new TH1D("", ";x [#mum]; x angle",
			pitchX, 0, pitchX);
	h_angleProfileX_rescuts_hits = new TH1D("", ";x [#mum]; x angle",
			pitchX, 0, pitchX);
	h_angleProfileY_rescuts_avg = new TH1D("", ";y [#mum]; y angle",
			pitchY, 0, pitchY);
	h_angleProfileY_rescuts_hits = new TH1D("", ";y [#mum]; y angle",
			pitchY, 0, pitchY);
}

/// do analysis of event for HIT ANALYSIS
void Botho::hit_analysis(const TbConfig &config, const Event &event) {
	// raw hits
	for(vector<PllHit *>::const_iterator it=event.rawHits.begin();
			it != event.rawHits.end(); it++) {
		int col = (*it)->col;
		int row = (*it)->row;
		h_rawHits->Fill(col, row);
	}
	
	// hit profiles
	h_profileX->Fill(event.trackX);
	h_profileY->Fill(event.trackY);
//	h_profileX_pixel->Fill(tbutils::getFoldedX(event));
	h_profileX_pixel->Fill(tbutils::getPixelX(event));	// (2014-10-20)
	h_profileY_pixel->Fill(tbutils::getPixelY(event));
	
	// hit map
	h_hitmap_raw_pos->Fill(event.trackX, event.trackY);
	
	// beam profiles
	h_beamProfileFEI3->Fill(event.trackX, event.trackY);
	h_beamProfileFEI4->Fill(event.trackX, event.trackY);
//	h_beamProfile2->Fill(event.trackX, event.trackY);
	
	// angle profile
	double dxdz = event.dxdz;
	double dydz = event.dydz;
	h_angleProfile->Fill(event.trackX, event.trackY,
			sqrt(dxdz*dxdz + dydz*dydz));
	h_angleProfile_numEvents->Fill(event.trackX, event.trackY);
	

	/// CUT 1:  EVENT.FTRACK == EVENT::KGOOD ///
	
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {  // if not event.fTrack == event::kGood
		return;                         // leave event
	}
		
	// hit multiplicity histogram
	h_hit_multiplicity->Fill(event.hits.size());
	
	// hit iden histogram
	for(vector<PllHit*>::const_iterator it = event.hits.begin();
			it != event.hits.end(); it++) {
		h_hit_iden->Fill((*it)->iden);
	}
			
	// hit maps
	for(vector<PllHit*>::const_iterator it = event.hits.begin();
			it != event.hits.end(); it++) {
		h_hitmap_all->Fill((*it)->col, (*it)->row);
	}
	// hit map all check
	for(int i = 0; i < event.hits.size(); i++) {
		h_hitmap_all_check->Fill(event.hits.at(i)->col, event.hits.at(i)->row);
	}
	h_hitmap_all_pos->Fill(event.trackX, event.trackY);
	
	// goodtrack profiles
	// fTrackRegion == event::kGood for normal FE-I3 modules
	// fTrackRegion is bad when either fTrackMaskedRegion or fTrackCentralRegion is bad
	// for FE-I4 modules differing conditions based on beam area
//	if(event.fTrackRegion == event::kGood) {  // commented out to make change 2014-10-24
	if((this->*track_geometry_and_mask)(event) == event::kGood) {
		h_profileX_pixel_goodtrackcuts->Fill(tbutils::getPixelX(event));
		h_profileY_pixel_goodtrackcuts->Fill(tbutils::getPixelY(event));
		h_profileX_module_goodtrackcuts->Fill(event.trackX);
		h_profileY_module_goodtrackcuts->Fill(event.trackY);
	}
	

	/// CUT 2: GET MATCHED CLUSTER ///
	
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
	                       // of matched cluster in vector of clusters of event
	if (matched == -1) {  // in case there is no matched cluster
		return;	// leave event
	}
		
	// get reference to the matched cluster
	const vector<PllHit*>& matched_cluster = event.clusters.at(matched);
	
	// hit multiplicity histogram matched
	h_hit_multiplicity_matched->Fill(event.hits.size());
	
	// matched hitmap positions
	h_hitmap_matched_pos->Fill(event.trackX, event.trackY);	
	
	
	/// CUT 3:  EVENT.FTRACKREGION == EVENT::KGOOD ///
	
	// Track in good region? (central and unmasked) (see event.h)
	//if (event.fTrackRegion != event::kGood) {
	// changed 2014-10-24
     if((this->*track_geometry_and_mask)(event) != event::kGood) {
		return;
	}
	
	// residuals
	double posX = cluster::getX(matched_cluster, event, cluster::kUnweighted);
	double posY = cluster::getY(matched_cluster, event, cluster::kUnweighted);
	  
	if (!(posX<0 || posY<0)){
		h_resX->Fill(event.trackX - posX);
		h_resY->Fill(event.trackY - posY);
	}
	
	int spanCols = cluster::spannedCols(matched_cluster);
	int spanRows = cluster::spannedRows(matched_cluster);
	
	// residuals for specific cluster sizes
	// CS1
	if(spanCols == 1)
		h_resX_CS1_unweighted->Fill(event.trackX - posX);
	if(spanRows == 1)
		h_resY_CS1_unweighted->Fill(event.trackY - posY);
	// CS2
	if(spanCols == 2)
		h_resX_CS2_unweighted->Fill(event.trackX - posX);
	if(spanRows == 2)
		h_resY_CS2_unweighted->Fill(event.trackY - posY);
/*	CS1
	if(matched_cluster.size() == 1) {
		h_resX_CS1_unweighted->Fill(event.trackX - posX);
		h_resY_CS1_unweighted->Fill(event.trackY - posY);
	}
	// CS2
	if(matched_cluster.size() == 2) {
		h_resX_CS2_unweighted->Fill(event.trackX - posX);
		h_resY_CS2_unweighted->Fill(event.trackY - posY);
	}
*/
	
	// fill MaxTotCell residual for all clustersizes
	posX = cluster::getX(matched_cluster, event, cluster::kMaxCell);
	posY = cluster::getY(matched_cluster, event, cluster::kMaxCell);
	
	h_resX_CSall_MaxTotCell->Fill(event.trackX - posX);
	h_resY_CSall_MaxTotCell->Fill(event.trackY - posY);

	if(spanRows == 1)
		h_resX_CSall_CXonly_MaxTotCell->Fill(event.trackX - posX);
	if(spanCols == 1)
		h_resY_CSall_CYonly_MaxTotCell->Fill(event.trackY - posY);
	
	// maxTOTcell residuals for restricted clusters
	if(spanRows <= 2 && spanCols <= 2) {
		h_resX_CS2x2max_MaxTotCell->Fill(event.trackX-posX);
		h_resY_CS2x2max_MaxTotCell->Fill(event.trackY-posY);
	}
	
	// pixel hit profiles for smoothing residuals
	h_profileX_pixel_rescuts->Fill(tbutils::getPixelX(event));
	h_profileY_pixel_rescuts->Fill(tbutils::getPixelY(event));
	
	// module hit profiles
	h_profileX_module_rescuts->Fill(event.trackX);
	h_profileY_module_rescuts->Fill(event.trackY);
	
	// angle pixel profile
//	h_angleProfileX_rescuts_avg->Fill(tbutils::getFoldedX
}

/// save histograms and free memory
void Botho::finalize_hit_analysis(const TbConfig &config) {
	// idens of hits (just for checking, whether all the hits in one analysis
	// actually have the same DUT iden (it is the DUT number specified in driver.cc))
	saveToFile(config, this->name, "hit_analysis", "hit_iden", h_hit_iden);
	
	// hit multiplicities
	saveToFile(config, this->name, "hit_analysis", "hit_multiplicity",
			h_hit_multiplicity);
	saveToFile(config, this->name, "hit_analysis", "hit_multiplicity_matched",
			h_hit_multiplicity_matched);
	
	// hitmaps
	saveToFile(config, this->name, "hit_analysis", "rawHits", h_rawHits);
	saveToFile(config, this->name, "hit_analysis", "hitmap_all", h_hitmap_all);
	saveToFile(config, this->name, "hit_analysis", "hitmap_all_check",
			h_hitmap_all_check);
	saveToFile(config, this->name, "hit_analysis", "hitmap_raw_pos", h_hitmap_raw_pos);
	saveToFile(config, this->name, "hit_analysis", "hitmap_all_pos", h_hitmap_all_pos);
	saveToFile(config, this->name, "hit_analysis", "hitmap_matched_pos",
			h_hitmap_matched_pos);
	
	// hit profiles
	saveToFile(config, this->name, "hit_analysis", "profileX", h_profileX);
	saveToFile(config, this->name, "hit_analysis", "profileY", h_profileY);
	saveToFile(config, this->name, "hit_analysis", "profileX_pixel",
			h_profileX_pixel);
	saveToFile(config, this->name, "hit_analysis", "profileY_pixel",
			h_profileY_pixel);
	
	// beam profile (different plots could be done with rebinning
	saveToFile(config, this->name, "hit_analysis", "beamProfileFEI3", h_beamProfileFEI3);
	saveToFile(config, this->name, "hit_analysis", "beamProfileFEI4", h_beamProfileFEI4);
//	saveToFile(config, this->name, "hit_analysis", "beamProfile2", h_beamProfile2);
	
	// angle profile
	h_angleProfile->Divide(h_angleProfile_numEvents);	// average angle profile
	saveToFile(config, this->name, "hit_analysis", "angleProfile", h_angleProfile);
 	
 	// residuals
 	saveToFile(config, this->name, "hit_analysis", "resX", h_resX);
 	saveToFile(config, this->name, "hit_analysis", "resY", h_resY);
 	
 	// residuals for specific cluster sizes
 	saveToFile(config, this->name, "hit_analysis", "resX_CS1_unweighted", h_resX_CS1_unweighted);
 	saveToFile(config, this->name, "hit_analysis", "resY_CS1_unweighted", h_resY_CS1_unweighted);
 	saveToFile(config, this->name, "hit_analysis", "resX_CS2_unweighted", h_resX_CS2_unweighted);
 	saveToFile(config, this->name, "hit_analysis", "resY_CS2_unweighted", h_resY_CS2_unweighted);
 	
 	saveToFile(config, this->name, "hit_analysis", "resX_CSall_MaxTotCell", h_resX_CSall_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resY_CSall_MaxTotCell", h_resY_CSall_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resX_CSall_CXonly_MaxTotCell", h_resX_CSall_CXonly_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resY_CSall_CYonly_MaxTotCell", h_resY_CSall_CYonly_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resX_CS2x2max_MaxTotCell", h_resX_CS2x2max_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resY_CS2x2max_MaxTotCell", h_resY_CS2x2max_MaxTotCell);
 	
 	// pixel hit profiles with same cuts as residuals
	saveToFile(config, this->name, "hit_analysis", "profileX_pixel_rescuts", h_profileX_pixel_rescuts);
	saveToFile(config, this->name, "hit_analysis", "profileY_pixel_rescuts", h_profileY_pixel_rescuts);
	
	saveToFile(config, this->name, "hit_analysis", "profileX_pixel_goodtrackcuts", h_profileX_pixel_goodtrackcuts);
	saveToFile(config, this->name, "hit_analysis", "profileY_pixel_goodtrackcuts", h_profileY_pixel_goodtrackcuts);
 	
 	// make factoring function from RAW profile BEFORE CUTS (bad english! I have a headake)
 	delete h_profileX_pixel_rescuts;
 	delete h_profileY_pixel_rescuts;
 	h_profileX_pixel_rescuts = h_profileX_pixel_goodtrackcuts;
 	h_profileY_pixel_rescuts = h_profileY_pixel_goodtrackcuts;
 	double max_valX = h_profileX_pixel_rescuts->GetMaximum();
 	double max_valY = h_profileY_pixel_rescuts->GetMaximum();
	h_profileX_pixel_rescuts->Sumw2();
 	h_profileX_pixel_rescuts->Scale(1.0/max_valX);
	h_profileY_pixel_rescuts->Sumw2();
 	h_profileY_pixel_rescuts->Scale(1.0/max_valY);
 	
 	TH1D *h_profileX_continuous = new TH1D("", ";x [#mum];factor", 4 * pitchX,
 			-2 * pitchX, 2*pitchX);
 	for(int i=0; i<=h_profileX_continuous->GetNbinsX(); i++) {
	 	h_profileX_continuous->SetBinContent(i+1,
	 			h_profileX_pixel_rescuts->GetBinContent(((i+(int)(pitchX/2.0))%((int)pitchX))+1));
	}
	h_resX_CS1_unweighted->Divide(h_profileX_continuous);
	h_resX_CSall_MaxTotCell->Divide(h_profileX_continuous);
	h_resX_CS2x2max_MaxTotCell->Divide(h_profileX_continuous);
	h_resX_CSall_CXonly_MaxTotCell->Divide(h_profileX_continuous);
	h_resX_CS1_unweighted->Sumw2();
	h_resX_CS1_unweighted->Scale(1.0/h_resX_CS1_unweighted->GetMaximum());
	h_resX_CSall_MaxTotCell->Sumw2();
	h_resX_CSall_MaxTotCell->Scale(1.0/h_resX_CSall_MaxTotCell->GetMaximum());
	h_resX_CS2x2max_MaxTotCell->Sumw2();
	h_resX_CS2x2max_MaxTotCell->Scale(1.0/h_resX_CS2x2max_MaxTotCell->GetMaximum());
	h_resX_CSall_CXonly_MaxTotCell->Sumw2();
	h_resX_CSall_CXonly_MaxTotCell->Scale(1.0/h_resX_CSall_CXonly_MaxTotCell->GetMaximum());
	saveToFile(config, this->name, "hit_analysis", "profileX_continuous_CS1", h_profileX_continuous);
	
	for(int i=0; i<=h_profileX_continuous->GetNbinsX(); i++) {
	 	h_profileX_continuous->SetBinContent(i+1,
	 			h_profileX_pixel_rescuts->GetBinContent((i%((int)pitchX))+1));
	}
	h_resX_CS2_unweighted->Divide(h_profileX_continuous);
	h_resX_CS2_unweighted->Sumw2();
	h_resX_CS2_unweighted->Scale(1.0/h_resX_CS2_unweighted->GetMaximum());
	saveToFile(config, this->name, "hit_analysis", "profileX_continuous_CS2", h_profileX_continuous);
	
 	TH1D *h_profileY_continuous = new TH1D("", ";x [#mum];factor", 4 * pitchY,
 			-2 * pitchY, 2*pitchY);
 	for(int i=0; i<=h_profileY_continuous->GetNbinsX(); i++) {
	 	h_profileY_continuous->SetBinContent(i+1,
	 			h_profileY_pixel_rescuts->GetBinContent(((i+(int)(pitchY/2.0))%((int)pitchY))+1));
	}
	h_resY_CS1_unweighted->Divide(h_profileY_continuous);
	h_resY_CSall_MaxTotCell->Divide(h_profileY_continuous);
	h_resY_CS2x2max_MaxTotCell->Divide(h_profileY_continuous);
	h_resY_CSall_CYonly_MaxTotCell->Divide(h_profileY_continuous);
	h_resY_CS1_unweighted->Sumw2();
	h_resY_CS1_unweighted->Scale(1.0/h_resY_CS1_unweighted->GetMaximum());
	h_resY_CSall_MaxTotCell->Sumw2();
	h_resY_CSall_MaxTotCell->Scale(1.0/h_resY_CSall_MaxTotCell->GetMaximum());
	h_resY_CS2x2max_MaxTotCell->Sumw2();
	h_resY_CS2x2max_MaxTotCell->Scale(1.0/h_resY_CS2x2max_MaxTotCell->GetMaximum());
	h_resY_CSall_CYonly_MaxTotCell->Sumw2();
	h_resY_CSall_CYonly_MaxTotCell->Scale(1.0/h_resY_CSall_CYonly_MaxTotCell->GetMaximum());
	saveToFile(config, this->name, "hit_analysis", "profileY_continuous_CS1", h_profileY_continuous);
	
	for(int i=0; i<=h_profileY_continuous->GetNbinsX(); i++) {
	 	h_profileY_continuous->SetBinContent(i+1,
	 			h_profileY_pixel_rescuts->GetBinContent((i%((int)pitchY))+1));
	}
	h_resY_CS2_unweighted->Divide(h_profileY_continuous);
	h_resY_CS2_unweighted->Sumw2();
	h_resY_CS2_unweighted->Scale(1.0/h_resY_CS2_unweighted->GetMaximum());
	saveToFile(config, this->name, "hit_analysis", "profileY_continuous_CS2", h_profileY_continuous);
	
	saveToFile(config, this->name, "hit_analysis", "resX_CS1_unweighted_normalized", h_resX_CS1_unweighted);
 	saveToFile(config, this->name, "hit_analysis", "resY_CS1_unweighted_normalized", h_resY_CS1_unweighted);
 	saveToFile(config, this->name, "hit_analysis", "resX_CS2_unweighted_normalized", h_resX_CS2_unweighted);
 	saveToFile(config, this->name, "hit_analysis", "resY_CS2_unweighted_normalized", h_resY_CS2_unweighted);
 	
 	saveToFile(config, this->name, "hit_analysis", "resX_CSall_MaxTotCell_normalized", h_resX_CSall_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resY_CSall_MaxTotCell_normalized", h_resY_CSall_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resX_CSall_CXonly_MaxTotCell_normalized", h_resX_CSall_CXonly_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resY_CSall_CYonly_MaxTotCell_normalized", h_resY_CSall_CYonly_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resX_CS2x2max_MaxTotCell_normalized", h_resX_CS2x2max_MaxTotCell);
 	saveToFile(config, this->name, "hit_analysis", "resY_CS2x2max_MaxTotCell_normalized", h_resY_CS2x2max_MaxTotCell);
 	
 	// module profile with same cuts as residuals
	saveToFile(config, this->name, "hit_analysis", "profileX_module_rescuts", h_profileX_module_rescuts);
	saveToFile(config, this->name, "hit_analysis", "profileY_module_rescuts", h_profileY_module_rescuts);
	// module profile with track kTrack.good und geometric cut used for residual correction
	saveToFile(config, this->name, "hit_analysis", "profileX_module_goodtrackcuts", h_profileX_module_goodtrackcuts);
	saveToFile(config, this->name, "hit_analysis", "profileY_module_goodtrackcuts", h_profileY_module_goodtrackcuts);
	 	
// 	for(int i=0; i<h_invProfile->GetNbinsX(); i++)
// 		h_invProfile->SetBinContent(i+1, max_val); 
 	
 	// free space of hit analysis histograms
 	delete h_hit_iden;
	
	delete h_hit_multiplicity;
	delete h_hit_multiplicity_matched;

	delete h_rawHits;	
	delete h_hitmap_all;
	delete h_hitmap_all_check;
	delete h_hitmap_raw_pos;
	delete h_hitmap_all_pos;
	
	delete h_profileX;
	delete h_profileY;
	delete h_profileX_pixel;
	delete h_profileY_pixel;
	
	delete h_beamProfileFEI3;
	delete h_beamProfileFEI4;
//	delete h_beamProfile2;
	
	delete h_angleProfile;
	delete h_angleProfile_numEvents;
	
	// residuals
	delete h_resX;
	delete h_resY;
	
	delete h_resX_CS1_unweighted;
	delete h_resY_CS1_unweighted;
	delete h_resX_CS2_unweighted;
	delete h_resY_CS2_unweighted;
	
	delete h_resX_CSall_MaxTotCell;
	delete h_resY_CSall_MaxTotCell;
	delete h_resX_CSall_CXonly_MaxTotCell;
	delete h_resY_CSall_CYonly_MaxTotCell;
	delete h_resX_CS2x2max_MaxTotCell;
	delete h_resY_CS2x2max_MaxTotCell;
	
//	delete h_profileX_pixel_rescuts;
//	delete h_profileY_pixel_rescuts;

	delete h_profileX_pixel_goodtrackcuts;
	delete h_profileY_pixel_goodtrackcuts;
	
	delete h_profileX_continuous;
	delete h_profileY_continuous;
	
	delete h_profileX_module_rescuts;
	delete h_profileY_module_rescuts;
	
	delete h_profileX_module_goodtrackcuts;
	delete h_profileY_module_goodtrackcuts;
}

////////////////////////////////////////////////////////////////////////////////
/// CLUSTER ANALYSIS        ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// initialize histograms etc
void Botho::init_cluster_analysis(TbConfig &config) {
	// copy DUT information from the config object to members of this object
	const DUT* dut = config.getDut(this->iden);	
	// get pitches
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();
	// get column and row numbers
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();
	
	// cluster size histograms
	h_clusize_all = new TH1D("", ";cluster size;entries",
			n_rows, +0.5, n_rows+0.5);
	h_clusize_matched = new TH1D("", ";cluster size;entries",
			n_rows, +0.5, n_rows+0.5);
	h_clusize_unmatched = new TH1D("", ";cluster size;entries",
			n_rows, +0.5, n_rows+0.5);
	h_clusize_matched_good_region = new TH1D("", ";cluster size;entries",
			n_rows, +0.5, n_rows+0.5);
	// cluster size histograms in directions
	h_clusize_x = new TH1D("", ";cluster size;entries",
			n_cols, +0.5, n_cols+0.5);
	h_clusize_x_CSx = new TH1D("", ";cluster size;entries",
			n_cols, +0.5, n_cols+0.5);
	h_clusize_y = new TH1D("", ";cluster size;entries",
			n_rows, +0.5, n_rows+0.5);
	h_clusize_y_CSy = new TH1D("", ";cluster size;entries",
			n_rows, +0.5, n_rows+0.5);
	// cluster multiplicity histograms
	h_cluster_multiplicity = new TH1D(
			"", ";number of clusters;events", 10, 0, 10);
	h_cluster_multiplicity_matched = new TH1D(
			"", ";number of clusters;events",10, 0, 10);		
	// cluster track maps
	for(int i=0; i<num_of_cs_track_histos; i++)
		h_track_cs[i] = new TH2D("", ";x [#mum];y [#mum]",
			(n_cols+4)*8, -2*pitchX, (n_cols+2)*pitchX,
			(n_rows+32) , -16*pitchY, (n_rows+16)*pitchY);
	// cluster track maps projected into single pixel
	for(int i=0; i<num_of_cs_track_histos; i++)
		h_pixel_track_cs[i] = new TH2D("", ";x [#mum];y [#mum]",
			80, 0, pitchX,
			10, 0, pitchY);
	h_pixel_track_cs_all = new TH2D("", "", 80, 0, pitchX,
			10, 0, pitchY);
	// cluster track maps projected into single pixel 1x1um binwidth
	for(int i=0; i<num_of_cs_track_histos; i++)
		h_pixel_track_1x1um_cs[i] = new TH2D("", ";x [#mum];y [#mum]",
			pitchX, 0, pitchX,
			pitchY, 0, pitchY);
	h_pixel_track_1x1um_cs_all = new TH2D("", "", pitchX, 0, pitchX,
			pitchY, 0, pitchY);
			
			
	// two times pixel pitch plots
	// for all cluster sizes
	h_pixel_track_cs_all_twopix = new TH2D("", "",
			TMath::Nint(pitchX)*2, 0, pitchX*2.0,
			TMath::Nint(pitchY)*2, 0, pitchY*2.0);
	// for cluster sizes 1 to num_of_cs_track_histos-1
	for(int i=0; i<num_of_cs_track_histos; i++) {
		h_pixel_track_cs_twopix[i] = new TH2D("", ";x [#mum];y [#mum]",
			TMath::Nint(pitchX)*2, 0, pitchX*2.0,
			TMath::Nint(pitchY)*2, 0, pitchY*2.0);
	}
	
	// high cluster size maps
	h_cs09_map = new TH2D("", ";Column;Row",
			n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_cs10_map = new TH2D("", ";Column;Row",
			n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_cs11_map = new TH2D("", ";Column;Row",
			n_cols, 0, n_cols,
			n_rows, 0, n_rows);
}
	
/// do analysis of event of CLUSTER ANALYSIS
void Botho::cluster_analysis(const TbConfig &config, const Event &event) {
	/// CUT 1:  EVENT.FTRACK == EVENT::KGOOD ///
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {  // if not event.fTrack == event::kGood
		return;                         // leave event
	}
	
	// cluster multiplicity histogram
	h_cluster_multiplicity->Fill(event.clusters.size());	// increase entry of
                                        // corresponding number of clusters by 1
	
	// cluster size all histogram
	for(vector<event::cluster_t>::const_iterator it = event.clusters.begin(); 
	          it != event.clusters.end(); it++)	{ // loop through all clusters 
	                                               // in event
		h_clusize_all->Fill(it->size());
	}
	// unmatched cluster size histogram
	for(vector<event::cluster_t>::const_iterator it = event.clusters.begin();
	     it != event.clusters.end(); it++) {	// loop over all clusters in event
		h_clusize_unmatched->Fill(it->size());	// unmatched histogram (matched
		                                   // ones will be removed after loop)
	}
	
	/// CUT 2: GET MATCHED CLUSTER ///
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
					// of matched cluster in vector of clusters of event
	if (matched == -1) {  // in case there is no matched cluster
		return;	// leave event
	}
	// get reference to the matched cluster
	const vector<PllHit*>& matched_cluster = event.clusters.at(matched);
	
	// matched cluster multiplicity histogram
	h_cluster_multiplicity_matched->Fill(event.clusters.size());
	
	// fill histogram cluster size matched
	h_clusize_matched->Fill(matched_cluster.size());
	// decrease unmatched histogram by value of matched cluster
	h_clusize_unmatched->Fill(matched_cluster.size(),-1.0);

	// high cluster size maps
	for(vector<PllHit*>::const_iterator it = matched_cluster.begin();
			it != matched_cluster.end(); it++) {
		if(matched_cluster.size()==9)
			h_cs09_map->Fill((*it)->col, (*it)->row);
		if(matched_cluster.size()==10)
			h_cs10_map->Fill((*it)->col, (*it)->row);
		if(matched_cluster.size()==11)
			h_cs11_map->Fill((*it)->col, (*it)->row);
	}

/*	
	/// CUT 3:  EVENT.FTRACKREGION == EVENT::KGOOD ///
	// Track in good region? (central and unmasked) (see event.h)
	if (event.fTrackRegion != event::kGood) {
		return;
	}
*/
	/// new CUT 3: introduced (2014-07-14) based on actual design of modules
	// tested by Botho
	// first check if track is inside masked pixel
	if((this->*track_check_central_region)(event) != event::kGood) {
		return;
	}
/*	if(event.fTrackMaskedRegion != event::kGood)
		return;
	if(isFEI3) {	// for FE-I3 cut top 7 rows and outer 400 um of first and last column
		if(event.trackX < 0 || event.trackX > 6800 ||
		   event.trackY > 7600) // event.trackY > (160-7)*50-25 (-25 safety)
			return;
	} else if(isFEI4) {
		if(event.fTrackRegion != event::kGood)
			return;
	} else {
		if(event.fTrackRegion != event::kGood)
			return;
	}
*/
	// module cluster track maps
	if(matched_cluster.size() < num_of_cs_track_histos) {	// check if cluster size smaller than constant
		if(matched_cluster.size() > 0) {	// check whether cluster size positive
			// fill cluster size histograms
			h_track_cs[matched_cluster.size()-1]->Fill(event.trackX, event.trackY);
		}
		else {	// if cluster size not greater than 0
			// debug message, there are no cluster size 0 expected
			TBALOG(kDEBUG2) <<
					"Matched event with cluster size"
					"less or equal to 0  does not make sense" << endl;
		}
	}	// fill overflow histogram for cluster size >= num_of_cs_track_histos
	else if(matched_cluster.size() >= num_of_cs_track_histos) {
		h_track_cs[num_of_cs_track_histos-1]->Fill(event.trackX, event.trackY);
	}
	
	// fill two pixel cluster size track maps
	h_pixel_track_cs_all_twopix->Fill(
			tbutils::getTrackTwoPixelEdgeDistance(event, "x"),
			tbutils::getTrackTwoPixelEdgeDistance(event, "y"));
	if(matched_cluster.size() < num_of_cs_track_histos) {	// check if cluster size smaller than constant
		if(matched_cluster.size() > 0) {	// check whether cluster size positive
			// fill cluster size histograms
			h_pixel_track_cs_twopix[matched_cluster.size()-1]->Fill(
					tbutils::getTrackTwoPixelEdgeDistance(event, "x"),
					tbutils::getTrackTwoPixelEdgeDistance(event, "y"));
		}
		else {	// if cluster size not greater than 0
			// debug message, there are no cluster size 0 expected
			TBALOG(kDEBUG2) <<
					"Matched event with cluster size"
					"less or equal to 0  does not make sense" << endl;
		}
	}	// fill overflow histogram for cluster size >= num_of_cs_track_histos
	else if(matched_cluster.size() >= num_of_cs_track_histos) {
		h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Fill(
				tbutils::getTrackTwoPixelEdgeDistance(event, "x"),
				tbutils::getTrackTwoPixelEdgeDistance(event, "y"));
	}
	
	// cluster track maps projected into single pixel
	double trackX = event.trackX;
	double trackY = event.trackY;
	trackX = trackX-pitchX/2.0;	// shift track pitchX/2 to the left to have first regular (not edge) column start at (0,0)
	int even_odd_col =0;	// counter to see whether column is even or odd and has to be mirrored
	while(trackX > pitchX) {
		trackX -= pitchX;
		even_odd_col++;
	}
	if(even_odd_col%2) {
		trackX = pitchX-trackX;
	}
	trackY = trackY-pitchY/2.0;	// same shift for trackY
	while(trackY > pitchY) {
		trackY -= pitchY;
	}
	// fill all track map
	h_pixel_track_cs_all->Fill(trackX, trackY);
	h_pixel_track_1x1um_cs_all->Fill(trackX, trackY);
	if(matched_cluster.size() < num_of_cs_track_histos) {	// check if cluster size smaller than constant
		if(matched_cluster.size() > 0) {	// check whether cluster size positive
			// fill cluster size histograms
			h_pixel_track_cs[matched_cluster.size()-1]->Fill(trackX, trackY);
			h_pixel_track_1x1um_cs[matched_cluster.size()-1]->Fill(trackX, trackY);
		}
		else {	// if cluster size not greater than 0
			// debug message, there are no cluster size 0 expected
			TBALOG(kDEBUG2) <<
					"Matched event with cluster size"
					"less or equal to 0  does not make sense" << endl;
		}
	}	// fill overflow histogram for cluster size >= num_of_cs_track_histos
	else if(matched_cluster.size() >= num_of_cs_track_histos) {
		h_pixel_track_cs[num_of_cs_track_histos-1]->Fill(trackX, trackY);
		h_pixel_track_1x1um_cs[num_of_cs_track_histos-1]->Fill(trackX, trackY);
	}
	

////////////////////////////////////////
/// DISABLED   CUT 3A DISABLE 2014-10-25	
	////////////////////////////////////////////
	///// CUT 3A: STANDARD REGION CUT
	////////////////////////////////////////////
//	if(event.fTrackRegion != event::kGood)
//		return;

	// MOVED TO CUT 3A 2014-10-23	// REMOVED again 2014-10-25	
	// fill diagram cluster size matched and good region
	h_clusize_matched_good_region->Fill(event.clusters.at(matched).size());
		
	// modified from clusterchecker.cc
	int MinSizeX = 10000;
	int MaxSizeX = 0;
	int MinSizeY = 10000;
	int MaxSizeY = 0;

//	const vector<PllHit*>& = matched_cluster;
	
	for (vector<PllHit*>::const_iterator jj =
					matched_cluster.begin();
					jj != matched_cluster.end(); ++jj) {
		if ((*jj)->col > MaxSizeX) {
			MaxSizeX = (*jj)->col;
		}
		if ((*jj)->col < MinSizeX) {
			MinSizeX = (*jj)->col;
		}
		if ((*jj)->row > MaxSizeY) {
			MaxSizeY = (*jj)->row;
		}
		if ((*jj)->row < MinSizeY) {
			MinSizeY = (*jj)->row;
		}
	}
	
//	h_matchClusterSize->Fill(event.clusters.at(matched).size());
	int clusize_x = MaxSizeX - MinSizeX + 1;
	int clusize_y = MaxSizeY - MinSizeY + 1;
	h_clusize_x->Fill(clusize_x);
	h_clusize_y->Fill(clusize_y);
	
	if(matched_cluster.size() == clusize_x)
		h_clusize_x_CSx->Fill(clusize_x);
	if(matched_cluster.size() == clusize_y)
		h_clusize_y_CSy->Fill(clusize_y);
}

/// save histograms and free memory for CLUSTER ANALYSIS
void Botho::finalize_cluster_analysis(const TbConfig &config) {
	// cluster histograms
	saveToFile(config, this->name, "cluster_analysis", "clusize_all",
			h_clusize_all);
	saveToFile(config, this->name, "cluster_analysis", "clusize_matched",
			h_clusize_matched);
	saveToFile(config, this->name, "cluster_analysis", "clusize_unmatched",
			h_clusize_unmatched);
	saveToFile(config, this->name, "cluster_analysis",
			"clusize_matched_good_region", h_clusize_matched_good_region);
			
	saveToFile(config, this->name, "cluster_analysis",
			"clusize_x", h_clusize_x);
	saveToFile(config, this->name, "cluster_analysis",
			"clusize_x_CSx", h_clusize_x_CSx);
	saveToFile(config, this->name, "cluster_analysis",
			"clusize_y", h_clusize_y);
	saveToFile(config, this->name, "cluster_analysis",
			"clusize_y_CSy", h_clusize_y_CSy);
	
	// cluster multiplicities
	saveToFile(config, this->name, "cluster_analysis", "cluster_multiplicity",
	          h_cluster_multiplicity);
	saveToFile(config, this->name, "cluster_analysis",
			"cluster_multiplicity_matched", h_cluster_multiplicity_matched);
	          
	// cluster track maps
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "track_cs", i+1);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_track_cs[i]);
	}
	saveToFile(config, this->name, "cluster_analysis", "track_cs_overflow",
			h_track_cs[num_of_cs_track_histos-1]);
			
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_cs_all",
			h_pixel_track_cs_all);
			
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_1x1um_cs_all",
			h_pixel_track_1x1um_cs_all);
	
	// cluster track maps projected into single pixel
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "pixel_track_cs", i+1);
		
		//divide by h_track_cs_all to normalize
		int num_entries = h_pixel_track_cs[i]->GetEntries();
		h_pixel_track_cs[i]->Divide(h_pixel_track_cs_all);
		h_pixel_track_cs[i]->SetEntries(num_entries);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs[i]);
				
		sprintf(histoname,"%s%d", "pixel_track_1x1um_cs_raw", i+1);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_1x1um_cs[i]);
		sprintf(histoname,"%s%d", "pixel_track_1x1um_cs", i+1);
		//divide by h_track_cs_all to normalize
		num_entries = h_pixel_track_1x1um_cs[i]->GetEntries();
		h_pixel_track_1x1um_cs[i]->Divide(h_pixel_track_1x1um_cs_all);
		h_pixel_track_1x1um_cs[i]->SetEntries(num_entries);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_1x1um_cs[i]);
	}
	int num_entries = h_pixel_track_cs[num_of_cs_track_histos-1]->GetEntries();
	h_pixel_track_cs[num_of_cs_track_histos-1]->Divide(h_pixel_track_cs_all);
	h_pixel_track_cs[num_of_cs_track_histos-1]->SetEntries(num_entries);
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_cs_overflow",
			h_pixel_track_cs[num_of_cs_track_histos-1]);
	num_entries = h_pixel_track_1x1um_cs[num_of_cs_track_histos-1]->GetEntries();
	h_pixel_track_1x1um_cs[num_of_cs_track_histos-1]->Divide(h_pixel_track_1x1um_cs_all);
	h_pixel_track_1x1um_cs[num_of_cs_track_histos-1]->SetEntries(num_entries);
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_1x1um_cs_overflow",
			h_pixel_track_1x1um_cs[num_of_cs_track_histos-1]);
			
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		// rebin cluster track maps projected into single pixel
		h_pixel_track_cs[i]->Rebin2D(2,2);
		h_pixel_track_cs[i]->Scale(1.0/4.0);
		char histoname[300];
		sprintf(histoname,"%s%d%s", "pixel_track_cs", i+1, "_rebin");
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs[i]);
	}
	h_pixel_track_cs[num_of_cs_track_histos-1]->Rebin2D(2,2);
	saveToFile(config, this->name, "cluster_analysis",
			"pixel_track_cs_overflow_rebin",
			h_pixel_track_cs[num_of_cs_track_histos-1]);
	
	// cluster track maps projected into two pixels
	// rebin all hit map
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_1x1um_cs_all_twopix",
			h_pixel_track_cs_all_twopix);
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "pixel_track_1x1um_cs_twopix_raw", i+1);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs_twopix[i]);
		sprintf(histoname,"%s%d", "pixel_track_1x1um_cs_twopix", i+1);
		
		// divide by h_track_cs_all to normalize
		int num_entries = h_pixel_track_cs_twopix[i]->GetEntries();
		h_pixel_track_cs_twopix[i]->Divide(h_pixel_track_cs_all_twopix);
		h_pixel_track_cs_twopix[i]->SetEntries(num_entries);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs_twopix[i]);
	}		
	// overflow map
	// divide by h_track_cs_all to normalize
	num_entries = h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->GetEntries();
	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Divide(h_pixel_track_cs_all_twopix);
	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->SetEntries(num_entries);
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_1x1um_cs_twopix_overflow",
			h_pixel_track_cs_twopix[num_of_cs_track_histos-1]);
			
				
	// rebin everything		
	h_pixel_track_cs_all_twopix->Rebin2D(cs_pix_maps_res_x, cs_pix_maps_res_y);
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_cs_all_twopix",
			h_pixel_track_cs_all_twopix);
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "pixel_track_cs_twopix", i+1);
		
		// rebin map
		h_pixel_track_cs_twopix[i]->Rebin2D(cs_pix_maps_res_x, cs_pix_maps_res_y);
		// divide by h_track_cs_all to normalize
//		int num_entries = h_pixel_track_cs_twopix[i]->GetEntries();
//		h_pixel_track_cs_twopix[i]->Divide(h_pixel_track_cs_all_twopix);
//		h_pixel_track_cs_twopix[i]->SetEntries(num_entries);
		h_pixel_track_cs_twopix[i]->Scale(1.0/(cs_pix_maps_res_x*cs_pix_maps_res_y));
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs_twopix[i]);
	}
	// overflow map
	// rebin map
	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Rebin2D(cs_pix_maps_res_x, cs_pix_maps_res_y);
	// divide by h_track_cs_all to normalize
//	num_entries = h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->GetEntries();
//	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Divide(h_pixel_track_cs_all_twopix);
//	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->SetEntries(num_entries);
	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Scale(1.0/(cs_pix_maps_res_x*cs_pix_maps_res_y));
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_cs_twopix_overflow",
			h_pixel_track_cs_twopix[num_of_cs_track_histos-1]);
			
	// rebin once more
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d_%dx%dum2", "pixel_track_cs_twopix", i+1,
			cs_pix_maps_res_x*2, cs_pix_maps_res_y*2);
		
		// rebin map
		h_pixel_track_cs_twopix[i]->Rebin2D(2, 2);
		h_pixel_track_cs_twopix[i]->Scale(1.0/4.0);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs_twopix[i]);
	}
	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Rebin2D(2, 2);
	h_pixel_track_cs_twopix[num_of_cs_track_histos-1]->Scale(4.0);
	char histoname[300];
	sprintf(histoname,"%s_%dx%dum2", "pixel_track_cs_twopix_overflow",
			cs_pix_maps_res_x*2, cs_pix_maps_res_y*2);
	saveToFile(config, this->name, "cluster_analysis", histoname,
			h_pixel_track_cs_twopix[num_of_cs_track_histos-1]);
			
	// high cluster size cluster maps
	saveToFile(config, this->name, "cluster_analysis", "cs09_map", h_cs09_map);
	saveToFile(config, this->name, "cluster_analysis", "cs10_map", h_cs10_map);
	saveToFile(config, this->name, "cluster_analysis", "cs11_map", h_cs11_map);
			
	// delete cluster histograms
	delete h_clusize_all;
	delete h_clusize_matched;
	delete h_clusize_unmatched;
	delete h_clusize_matched_good_region;
	
	delete h_clusize_x;
	delete h_clusize_x_CSx;
	delete h_clusize_y;
	delete h_clusize_y_CSy;
	
	delete h_cluster_multiplicity;
	delete h_cluster_multiplicity_matched;
	
	for(int i=0; i < num_of_cs_track_histos; i++)
		delete h_track_cs[i];
	delete h_pixel_track_cs_all;
	delete h_pixel_track_1x1um_cs_all;
	for(int i=0; i < num_of_cs_track_histos; i++) {
		delete h_pixel_track_cs[i];
		delete h_pixel_track_1x1um_cs[i];
	}
	delete h_pixel_track_cs_all_twopix;
	for(int i=0; i < num_of_cs_track_histos; i++) {
		delete h_pixel_track_cs_twopix[i];
	}
		
	delete h_cs09_map;
	delete h_cs10_map;
	delete h_cs11_map;
}


////////////////////////////////////////////////////////////////////////////////
/// TOT - CHARGE ANALYSIS   ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// initialize tot analysis
void Botho::init_tot_analysis(TbConfig &config) {
	// copy DUT information from the config object to members of this object
	const DUT* dut = config.getDut(this->iden);
	
	// get pitches
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();
	
	// get column and row numbers
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();
	
	// number of tot bins
	int n_tot_bins; 
	if(n_cols==80 && n_rows==336) { // case of FE-I4
		n_tot_bins = 16;
	} 
	else if(n_cols==18 && n_rows==160) {	// case of FE-I3
		n_tot_bins = 256;
	} 
	else if(n_cols==40 && n_rows==672) {// case of FE-I4-25
	        n_tot_bins = 16;
	}
	else{
		TBALOG( kERROR ) << "unknown type of device"
				<< " -> cannot assign ToT bin number"
				<< std::endl << "exiting" << std::endl;
		std::exit(-1);
	}
	
	/// ToT histograms
	h_tot_all = new TH1D("", ";ToT;Entries", 3*n_tot_bins, 0, 3*n_tot_bins);
	h_tot_cs1 = new TH1D("", ";ToT;Entries", n_tot_bins, 0, n_tot_bins);
	h_tot_cs2 = new TH1D("", ";ToT;Entries", 2*n_tot_bins, 0, 2*n_tot_bins);
	
	/// ToT module map
	h_module_tot_avg = new TH2D("", ";x;y;ToT", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_module_tot_hits = new TH2D("", ";x;y;hits", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	// ToT module map cluster size 1 hits
	h_module_tot_avg_cs1 = new TH2D("", ";x;y;ToT", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_module_tot_hits_cs1 = new TH2D("", ";x;y;hits", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	// ToT module map cluster size 2 hits
	h_module_tot_avg_cs2 = new TH2D("", ";x;y;ToT", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_module_tot_hits_cs2 = new TH2D("", ";x;y;hits", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
		
			
	/// ToT pixel maps
	int xRes = 2;
	int yRes = 2;
	int res1 = 5;
	h_pixel_tot_avg = new TH2D("", ";x [#mum];y [#mum];ToT",
			pitchX/xRes, 0, pitchX, pitchY/yRes, 0, pitchY);
	h_pixel_tot_hits = new TH2D("", ";x [#mum];y [#mum];hits",
			pitchX/xRes, 0, pitchX, pitchY/yRes, 0, pitchY);
	h_pixel_tot_avg_alt = new TH2D("", ";x [#mum];y [#mum];ToT",
			pitchX/res1, 0, pitchX, pitchY/res1, 0, pitchY);
	h_pixel_tot_hits_alt = new TH2D("", ";x [#mum];y [#mum];hits",
			pitchX/res1, 0, pitchX, pitchY/res1, 0, pitchY);
			
			
	/// ToT pixel maps of certain CS
	for(int i=0; i<4; i++) {
		h_pixel_tot_avg_CS[i] = new TH2D("", ";x [#mum];y [#mum];ToT",
			pitchX/xRes, 0, pitchX, pitchY/yRes, 0, pitchY);
		h_pixel_tot_hits_CS[i] = new TH2D("", ";x [#mum];y [#mum];ToT",
			pitchX/xRes, 0, pitchX, pitchY/yRes, 0, pitchY);
	}
			
			
	/// ToT 3D scatter histogram
	res1 = 5;	// 5x5 um resolution
	h_pixel_tot_3d = new TH3D("", ";x [#mum];y [#mum];ToT",
			pitchX/res1, 0, pitchX, pitchY/res1, 0, pitchY,
			150, 0, 150);
			
	create_totmap_file =
			(bool) config.cmdLineExtras_argGetter(
					(string) "A_botho_CreateTotmapFile", (bool) false,
					(string) "switches on creating a totmap file");
					
	do_totmap_selection =
			(bool) config.cmdLineExtras_argGetter(
					(string) "A_botho_DoTotmapSelection", (bool) false,
					(string) "switches on doing totmap depending analysis");
					
	do_square_selection =
			(bool) config.cmdLineExtras_argGetter(
					(string) "A_botho_DoSquareSelection", (bool) false,
					(string) "switches on doing square area depending analysis");
	
	// if square area selction is activated
	if(do_square_selection) {
		// initialize tot histograms
		for(int i=0; i < 2; i++) {
			h_tot_all_square_selected[i] = new TH1D("", ";ToT;Entries",
					3*n_tot_bins, 0, 3*n_tot_bins);
			h_tot_cs1_square_selected[i] = new TH1D("", ";ToT;Entries",
					n_tot_bins, 0, n_tot_bins);
			h_tot_cs2_square_selected[i] = new TH1D("", ";ToT;Entries",
					2*n_tot_bins, 0, 2*n_tot_bins);
		}
		
		// initialize division maps
		for(int i=0; i < 2; i++) {
			h_module_totmap_square_selected[i] = new TH2D(*h_module_tot_avg);
		}
		// loop over all pixels
		for(int x=0; x < h_module_tot_avg->GetNbinsX(); x++) {
			for(int y=0; y < h_module_tot_avg->GetNbinsY(); y++) {
				if(in_square(ROI_square, x, y)) {
					h_module_totmap_square_selected[0]->Fill(x,y);
				} else {
					h_module_totmap_square_selected[1]->Fill(x,y);
				}
			}
		}
		
		
		// initialize pixel map selected
		for(int i=0; i < 2; i++) {
			h_pixel_tot_hits_square_selected[i]
					= new TH2D("", ";x [#mum];y [#mum];ToT",
					pitchX, 0, pitchX, pitchY, 0, pitchY);
			h_pixel_tot_avg_square_selected [i]
					= new TH2D("", ";x [#mum];y [#mum];ToT",
					pitchX, 0, pitchX, pitchY, 0, pitchY);
		}
	}
	
	
	// if totmap selection is activated
	if(do_totmap_selection) {
		// calculate of span of tot selection
		tot_sel_span = std::fabs(tot_selected_upperbound
				- tot_selected_lowerbound);
		// calculate size of one "bin"
		tot_sel_division = tot_sel_span/num_of_tot_selected_divisions;
		// initialize h_module_tot_selection as a copy of h_module_tot_avg
		h_module_tot_selection = new TH2D(*h_module_tot_avg);
		
		// initialize tot histograms
		for(int i=0; i < num_of_tot_selected_divisions+2; i++) {
			h_tot_all_selected[i] = new TH1D("", ";ToT;Entries",
					3*n_tot_bins, 0, 3*n_tot_bins);
			h_tot_cs1_selected[i] = new TH1D("", ";ToT;Entries",
					n_tot_bins, 0, n_tot_bins);
			h_tot_cs2_selected[i] = new TH1D("", ";ToT;Entries",
					2*n_tot_bins, 0, 2*n_tot_bins);
		}
		
		// initialize division maps (don't forget overflow and underflow)
		for(int i=0; i < num_of_tot_selected_divisions+2; i++) {
			h_module_totmap_selected[i] = new TH2D(*h_module_tot_avg);
		}
		
		// initialize pixel map selected
		for(int i=0; i < num_of_tot_selected_divisions+2; i++) {
			h_pixel_tot_hits_selected[i] = new TH2D(*h_pixel_tot_hits);
			h_pixel_tot_avg_selected [i] = new TH2D(*h_pixel_tot_avg);
		}
		
		// read in totmap selection file
		// the file should contain all average tot values as written by
		// botho::finalize function
		std::string line;
		// get file name from config function
		char filename[500];
		sprintf(filename, "%s", config.getOutStreamName(this->name, "totmap"));
		std::ifstream instream(filename);
		if (!instream.is_open()) {
			std::cerr << "The totmap file '" << filename
					<< "' could not be opened."
					<< std::endl;
			exit(-1);
		}
		// loop over all pixels
		for(int x=0; x < h_module_tot_avg->GetNbinsX(); x++) {
			for(int y=0; y < h_module_tot_avg->GetNbinsY(); y++) {
				if(std::getline(instream, line).eof()) {
					break;
				}
				// convert line to double value
				double value = std::strtod(line.c_str(), NULL);
				h_module_tot_selection->Fill(x, y, value);
				int i=0;
				while(i < num_of_tot_selected_divisions) {
				//for(int i=0;i < num_of_tot_selected_divisions; i++) {
					if((tot_selected_lowerbound
							+i*tot_sel_division <= value)
							&& (value < tot_selected_lowerbound
							+(i+1)*tot_sel_division)) {
						h_module_totmap_selected[i+1]->Fill(x,y);
						break;
					}
					i++;
				}
				if(i == num_of_tot_selected_divisions) {
					if(value < tot_selected_lowerbound) {
						h_module_totmap_selected[0]->Fill(x,y);
					}
					if(value >= tot_selected_upperbound) {
						h_module_totmap_selected[num_of_tot_selected_divisions+1]->Fill(x,y);
					}
				}
			}
		}
		instream.close();
	}	
}

/// do analysis of event for tot analysis
void Botho::tot_analysis(const TbConfig &config, const Event &event) {
	/// CUT 1:  EVENT.FTRACK == EVENT::KGOOD ///
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {  // if not event.fTrack == event::kGood
		return;                         // leave event
	}
	/// CUT 2: GET MATCHED CLUSTER ///
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
					// of matched cluster in vector of clusters of event
	if (matched == -1) {  // in case there is no matched cluster
		return;	// leave event
	}
	/// CUT 3:  EVENT.FTRACKREGION == EVENT::KGOOD ///
	// Track in good region? (central and unmasked) (see event.h)
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	
	// get reference to the matched cluster
	const vector<PllHit*>& matched_cluster = event.clusters.at(matched);
	
	// ToT histogram
	h_tot_all->Fill(cluster::getSumTot(matched_cluster));
	// ToT module histograms (unweighted algorithm for cluster position)
	h_module_tot_avg->Fill(
			cluster::getUnWeightedCol(matched_cluster),
			cluster::getUnWeightedRow(matched_cluster),
			cluster::getSumTot(matched_cluster));
	h_module_tot_hits->Fill(
			cluster::getUnWeightedCol(matched_cluster),
			cluster::getUnWeightedRow(matched_cluster));
	// ToT CS1
	if(matched_cluster.size() == 1) {
		h_tot_cs1->Fill(cluster::getSumTot(matched_cluster));
		h_module_tot_avg_cs1->Fill(
				cluster::getUnWeightedCol(matched_cluster),
				cluster::getUnWeightedRow(matched_cluster),
				cluster::getSumTot(matched_cluster));
		h_module_tot_hits_cs1->Fill(
				cluster::getUnWeightedCol(matched_cluster),
				cluster::getUnWeightedRow(matched_cluster));
	} else if(matched_cluster.size() == 2) {	// CS2
		h_tot_cs2->Fill(cluster::getSumTot(matched_cluster));
		h_module_tot_avg_cs2->Fill(
				cluster::getUnWeightedCol(matched_cluster),
				cluster::getUnWeightedRow(matched_cluster),
				cluster::getSumTot(matched_cluster));
		h_module_tot_hits_cs2->Fill(
				cluster::getUnWeightedCol(matched_cluster),
				cluster::getUnWeightedRow(matched_cluster));
	}
			
	// ToT pixel histograms (from stefano version of qEfficiency.cc)
	// with additional cuts
	if ((event.fTrack             		== event::kGood)	// cut 1
			&& (event.fClusters           == event::kGood)
			&& (event.hits.size()         != 0           )
			&& (event.fTrackCentralRegion == event::kGood)
	  		&& (event.fTrackRegion        == event::kGood)) {	// cut 3
		h_pixel_tot_avg->Fill(
				tbutils::getFoldedX(event), tbutils::getPixelY(event),
				cluster::getSumTot(matched_cluster));
		if(cluster::getSumTot(matched_cluster)>0){
			h_pixel_tot_hits->Fill(
					tbutils::getFoldedX(event), tbutils::getPixelY(event));
		}
		// fill ToT pixel maps with for CS
		for(int i=0; i<4; i++) {
			if(matched_cluster.size() == i+1) {
				h_pixel_tot_avg_CS[i]->Fill(
						tbutils::getFoldedX(event), tbutils::getPixelY(event),
						cluster::getSumTot(matched_cluster));
				h_pixel_tot_hits_CS[i]->Fill(
						tbutils::getFoldedX(event), tbutils::getPixelY(event));
			}
		}
		
		// fill selected pixel maps
		if(do_square_selection) {
			double value = cluster::getSumTot(matched_cluster);
			int i = 0;
			if(in_square(ROI_square,
					cluster::getUnWeightedCol(matched_cluster),
					cluster::getUnWeightedRow(matched_cluster))) {
				i = 0;
			} else {
				i = 1;
			}
			// ToT histograms
			h_tot_all_square_selected[i]->
					Fill(value);
			// ToT histograms CS1
			if(matched_cluster.size() == 1) {
				h_tot_cs1_square_selected[i]->
						Fill(value);
			} else if(matched_cluster.size() == 2) {	// CS2
				h_tot_cs2_square_selected[i]->
						Fill(value);
			}
			// ToT maps
			h_pixel_tot_hits_square_selected[i]->Fill(
					tbutils::getFoldedX(event),
					tbutils::getPixelY(event));
			h_pixel_tot_avg_square_selected [i]->Fill(
					tbutils::getFoldedX(event),
					tbutils::getPixelY(event),
					value);
		}
		
		// FIXME could probably be implemented more elegantly with the
		// 		already created h_module_tot_selected[] maps
		if(do_totmap_selection) {
			double sel_value = h_module_tot_selection->GetBinContent(
					h_module_tot_selection->FindBin(
					cluster::getUnWeightedCol(matched_cluster),
					cluster::getUnWeightedRow(matched_cluster)));
			double value = cluster::getSumTot(matched_cluster);
			int i=0;
			while(i < num_of_tot_selected_divisions) {
			//for(int i=0;i < num_of_tot_selected_divisions; i++) {
				if((tot_selected_lowerbound
						+i*tot_sel_division <= sel_value)
						&& (sel_value < tot_selected_lowerbound
						+(i+1)*tot_sel_division)) {
					// ToT histograms
					h_tot_all_selected[i+1]->
							Fill(value);
					// ToT histograms CS1
					if(matched_cluster.size() == 1) {
						h_tot_cs1_selected[i+1]->
								Fill(value);
					} else if(matched_cluster.size() == 2) {	// CS2
						h_tot_cs2_selected[i+1]->
								Fill(value);
					}
					// ToT maps
					h_pixel_tot_hits_selected[i+1]->Fill(
							tbutils::getFoldedX(event),
							tbutils::getPixelY(event));
					h_pixel_tot_avg_selected [i+1]->Fill(
							tbutils::getFoldedX(event),
							tbutils::getPixelY(event),
							value);
					break;
				}
				i++;
			}
			if(i == num_of_tot_selected_divisions) {
				if(sel_value < tot_selected_lowerbound) {
					// ToT histograms
					h_tot_all_selected[0]->
							Fill(value);
					// ToT histograms CS1
					if(matched_cluster.size() == 1) {
						h_tot_cs1_selected[0]->
								Fill(value);
					} else if(matched_cluster.size() == 2) {	// CS2
						h_tot_cs2_selected[0]->
								Fill(value);
					}
					// ToT maps
					h_pixel_tot_hits_selected[0]->Fill(
							tbutils::getFoldedX(event),
							tbutils::getPixelY(event));
					h_pixel_tot_avg_selected [0]->Fill(
							tbutils::getFoldedX(event),
							tbutils::getPixelY(event),
							value);
				}
				if(sel_value >= tot_selected_upperbound) {
					// ToT histograms
					h_tot_all_selected[num_of_tot_selected_divisions+1]->
							Fill(value);
					// ToT histograms CS1
					if(matched_cluster.size() == 1) {
						h_tot_cs1_selected[num_of_tot_selected_divisions+1]->
								Fill(value);
					} else if(matched_cluster.size() == 2) {	// CS2
						h_tot_cs2_selected[num_of_tot_selected_divisions+1]->
								Fill(value);
					}
					h_pixel_tot_hits_selected[num_of_tot_selected_divisions+1]->Fill(
							tbutils::getFoldedX(event),
							tbutils::getPixelY(event));
					h_pixel_tot_avg_selected [num_of_tot_selected_divisions+1]->Fill(
							tbutils::getFoldedX(event),
							tbutils::getPixelY(event),
							value);
				}
			}
		}
	}
	
	// alternative ToT pixel histograms by Botho
	h_pixel_tot_avg_alt->Fill(
			tbutils::getFoldedX(event), tbutils::getPixelY(event),
			cluster::getSumTot(matched_cluster));
	h_pixel_tot_hits_alt->Fill(
			tbutils::getFoldedX(event), tbutils::getPixelY(event));
			
	// 3D ToT scatter histogram
	h_pixel_tot_3d->Fill(tbutils::getFoldedX(event), tbutils::getPixelY(event),
			cluster::getSumTot(matched_cluster));
}

/// save histograms and free memory of tot analysis
void Botho::finalize_tot_analysis(const TbConfig &config) {
	// ToT histograms
	saveToFile(config, this->name, "tot_analysis", "tot_all",
			h_tot_all);
	saveToFile(config, this->name, "tot_analysis", "tot_cs1",
			h_tot_cs1);
	saveToFile(config, this->name, "tot_analysis", "tot_cs2",
			h_tot_cs2);

	// ToT module maps
	saveToFile(config, this->name, "tot_analysis", "module_tot_sum",
			h_module_tot_avg);
	h_module_tot_avg->Divide(h_module_tot_hits);	// calculate average ToT
	saveToFile(config, this->name, "tot_analysis", "module_tot_hits",
			h_module_tot_hits);
	saveToFile(config, this->name, "tot_analysis", "module_tot_avg",
			h_module_tot_avg);
	// ToT module maps cluster size 1 hits
	saveToFile(config, this->name, "tot_analysis", "module_tot_sum_cs1",
			h_module_tot_avg_cs1);
	h_module_tot_avg_cs1->Divide(h_module_tot_hits_cs1);	// calculate average ToT
	saveToFile(config, this->name, "tot_analysis", "module_tot_hits_cs1",
			h_module_tot_hits_cs1);
	saveToFile(config, this->name, "tot_analysis", "module_tot_avg_cs1",
			h_module_tot_avg_cs1);
	// ToT module maps cluster size 2 hits
	saveToFile(config, this->name, "tot_analysis", "module_tot_sum_cs2",
			h_module_tot_avg_cs2);
	h_module_tot_avg_cs2->Divide(h_module_tot_hits_cs2);	// calculate average ToT
	saveToFile(config, this->name, "tot_analysis", "module_tot_hits_cs2",
			h_module_tot_hits_cs2);
	saveToFile(config, this->name, "tot_analysis", "module_tot_avg_cs2",
			h_module_tot_avg_cs2);
	// write down TOT map to file in case the command line argument is given
	if(create_totmap_file) {
		std::ofstream outstream;
		// get filename from config function
		outstream.open(config.getOutStreamName(this->name, "totmap"));
		for(int x=1; x <= h_module_tot_avg->GetNbinsX(); x++) {
			for(int y=1; y <= h_module_tot_avg->GetNbinsY(); y++) {
				outstream << h_module_tot_avg->GetBinContent(x, y) << std::endl;
			}
		}
		outstream.close();
	}
	// results depending on square area
	if(do_square_selection) {
		// square selection maps
		saveToFile(config, this->name, "tot_analysis",
				"module_totmap_square_selected_0",
				h_module_totmap_square_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"module_totmap_square_selected_1",
				h_module_totmap_square_selected[1]);
		// histograms
		saveToFile(config, this->name, "tot_analysis",
				"tot_all_square_selected_0",
				h_tot_all_square_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_all_square_selected_1",
				h_tot_all_square_selected[1]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_cs1_square_selected_0",
				h_tot_cs1_square_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_cs1_square_selected_1",
				h_tot_cs1_square_selected[1]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_cs2_square_selected_0",
				h_tot_cs2_square_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_cs2_square_selected_1",
				h_tot_cs2_square_selected[1]);
		// pixel maps
		for(int i=0; i < 2; i++) {
			char name[200];
			// save pixel hit map
			sprintf(name,"pixel_tot_hits_square_selected_%i", i);
			saveToFile(config, this->name, "tot_analysis", name,
					h_pixel_tot_hits_square_selected[i]);
			// save pixel tot sum map
			sprintf(name,"pixel_tot_sum_square_selected_%i", i);
			saveToFile(config, this->name, "tot_analysis", name,
					h_pixel_tot_avg_square_selected[i]);
			// calculate and save pixel tot avg map
			h_pixel_tot_avg_square_selected[i]->
					Divide(h_pixel_tot_hits_square_selected[i]);
			sprintf(name,"pixel_tot_avg_square_selected_%i", i);
			saveToFile(config, this->name, "tot_analysis", name,
					h_pixel_tot_avg_square_selected[i]);
			// save more coarse rebinned versions
			h_pixel_tot_avg_square_selected[i]->Rebin2D(5,5);
			h_pixel_tot_avg_square_selected[i]->Scale(0.04); // scale 1/25
			char rebin_name[200];
			sprintf(rebin_name,"%s_rebin1",name);
			saveToFile(config, this->name, "tot_analysis", rebin_name,
					h_pixel_tot_avg_square_selected[i]);
			// second rebin
			h_pixel_tot_avg_square_selected[i]->Rebin2D(2,2);
			h_pixel_tot_avg_square_selected[i]->Scale(0.25); // scale 1/4
			sprintf(rebin_name,"%s_rebin2",name);
			saveToFile(config, this->name, "tot_analysis", rebin_name,
					h_pixel_tot_avg_square_selected[i]);
		}
	}
	
	// results depending on totmap value preselection
	if(do_totmap_selection) {
		saveToFile(config, this->name, "tot_analysis", "module_tot_selection",
				h_module_tot_selection);
		// ToT histograms
		for(int i=0; i < num_of_tot_selected_divisions; i++) {
			char basename[200];
			sprintf(basename,"tot_selected_%i_%i",
				TMath::Nint((tot_selected_lowerbound+i*tot_sel_division)*10),
				TMath::Nint((tot_selected_lowerbound+(i+1)*tot_sel_division)*10));
			char name[200];
			sprintf(name, "%s_all", basename);
			saveToFile(config, this->name, "tot_analysis", name,
					h_tot_all_selected[i+1]);
			sprintf(name, "%s_cs1", basename);
			saveToFile(config, this->name, "tot_analysis", name,
					h_tot_cs1_selected[i+1]);
			sprintf(name, "%s_cs2", basename);
			saveToFile(config, this->name, "tot_analysis", name,
					h_tot_cs2_selected[i+1]);
		}
		// underflow
		saveToFile(config, this->name, "tot_analysis",
				"tot_selected_underflow_all",
				h_tot_all_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_selected_underflow_cs1",
				h_tot_cs1_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_selected_underflow_cs2",
				h_tot_cs2_selected[0]);
		// overflow
		saveToFile(config, this->name, "tot_analysis",
				"tot_selected_overflow_all",
				h_tot_all_selected[num_of_tot_selected_divisions+1]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_selected_overflow_cs1",
				h_tot_cs1_selected[num_of_tot_selected_divisions+1]);
		saveToFile(config, this->name, "tot_analysis",
				"tot_selected_overflow_cs2",
				h_tot_cs2_selected[num_of_tot_selected_divisions+1]);
				
				
		for(int i=0; i < num_of_tot_selected_divisions; i++) {
			char name[200];
			sprintf(name,"module_tot_selected_%i_%i",
				TMath::Nint((tot_selected_lowerbound+i*tot_sel_division)*10),
				TMath::Nint((tot_selected_lowerbound+(i+1)*tot_sel_division)*10));
			saveToFile(config, this->name, "tot_analysis", name,
					h_module_totmap_selected[i+1]);
		}
		saveToFile(config, this->name, "tot_analysis",
				"module_tot_selected_underflow",
				h_module_totmap_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"module_tot_selected_overflow",
				h_module_totmap_selected[num_of_tot_selected_divisions+1]);
		for(int i=0; i < num_of_tot_selected_divisions; i++) {
			char name[200];
			// save pixel hit map
			sprintf(name,"pixel_tot_hits_selected_%i_%i",
				TMath::Nint((tot_selected_lowerbound+i*tot_sel_division)*10),
				TMath::Nint((tot_selected_lowerbound+(i+1)*tot_sel_division)*10));
			saveToFile(config, this->name, "tot_analysis", name,
					h_pixel_tot_hits_selected[i+1]);
			// save pixel tot sum map
			sprintf(name,"pixel_tot_sum_selected_%i_%i",
				TMath::Nint((tot_selected_lowerbound+i*tot_sel_division)*10),
				TMath::Nint((tot_selected_lowerbound+(i+1)*tot_sel_division)*10));
			saveToFile(config, this->name, "tot_analysis", name,
					h_pixel_tot_avg_selected[i+1]);
			// calculate and save pixel tot avg map
			h_pixel_tot_avg_selected[i+1]->Divide(h_pixel_tot_hits_selected[i+1]);
			sprintf(name,"pixel_tot_avg_selected_%i_%i",
				TMath::Nint((tot_selected_lowerbound+i*tot_sel_division)*10),
				TMath::Nint((tot_selected_lowerbound+(i+1)*tot_sel_division)*10));
			saveToFile(config, this->name, "tot_analysis", name,
					h_pixel_tot_avg_selected[i+1]);
			// save more coarse rebinned versions
			h_pixel_tot_avg_selected[i+1]->Rebin2D(2,2);
			h_pixel_tot_avg_selected[i+1]->Scale(0.25); // scale 1/4
			char rebin_name[200];
			sprintf(rebin_name,"%s_rebin1",name);
			saveToFile(config, this->name, "tot_analysis", rebin_name,
					h_pixel_tot_avg_selected[i+1]);
			// second rebin
			h_pixel_tot_avg_selected[i+1]->Rebin2D(2,2);
			h_pixel_tot_avg_selected[i+1]->Scale(0.25); // scale 1/4
			sprintf(rebin_name,"%s_rebin2",name);
			saveToFile(config, this->name, "tot_analysis", rebin_name,
					h_pixel_tot_avg_selected[i+1]);
		}
		// underflow pixel hits
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_hits_selected_underflow",
				h_pixel_tot_hits_selected[0]);
		// underflow pixel sum
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_sum_selected_underflow",
				h_pixel_tot_avg_selected[0]);
		// underflow pixel avg
		h_pixel_tot_avg_selected[0]->Divide(h_pixel_tot_hits_selected[0]);
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_avg_selected_underflow",
				h_pixel_tot_avg_selected[0]);
		// save more coarse rebinned versions
		h_pixel_tot_avg_selected[0]->Rebin2D(2,2);
		h_pixel_tot_avg_selected[0]->Scale(0.25); // scale 1/4
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_avg_selected_underflow_rebin1",
				h_pixel_tot_avg_selected[0]);
		// second rebin
		h_pixel_tot_avg_selected[0]->Rebin2D(2,2);
		h_pixel_tot_avg_selected[0]->Scale(0.25); // scale 1/4
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_avg_selected_underflow_rebin2",
				h_pixel_tot_avg_selected[0]);
		// overflow pixel hits
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_hits_selected_overflow",
				h_pixel_tot_hits_selected[num_of_tot_selected_divisions+1]);
		// overflow pixel sum
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_sum_selected_overflow",
				h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]);
		// overflow pixel avg
		h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]->
				Divide(h_pixel_tot_hits_selected[num_of_tot_selected_divisions+1]);
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_avg_selected_overflow",
				h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]);
		// save more coarse rebinned versions
		h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]->
				Rebin2D(2,2);
		h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]->
				Scale(0.25); // scale 1/4
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_sum_selected_overflow_rebin1",
				h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]);
		// second rebin
		h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]->
				Rebin2D(2,2);
		h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]->
				Scale(0.25); // scale 1/4
		saveToFile(config, this->name, "tot_analysis",
				"pixel_tot_sum_selected_overflow_rebin2",
				h_pixel_tot_avg_selected[num_of_tot_selected_divisions+1]);
	}
	
	// ToT pixel maps
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_sum",
			h_pixel_tot_avg);
	h_pixel_tot_avg->Divide(h_pixel_tot_hits);	// calculate average ToT
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_hits",
			h_pixel_tot_hits);
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_avg",
			h_pixel_tot_avg);
			
	// ToT pixel maps for CS
	for(int i; i<4; i++) {
		char name[100];
		sprintf(name,"pixel_tot_sum_CS%i", i+1);
		saveToFile(config, this->name, "tot_analysis", name,
				h_pixel_tot_avg_CS[i]);
		h_pixel_tot_avg_CS[i]->Divide(h_pixel_tot_hits_CS[i]);	// calculate average ToT
		sprintf(name,"pixel_tot_hits_CS%i", i+1);
		saveToFile(config, this->name, "tot_analysis", name,
				h_pixel_tot_hits_CS[i]);
		sprintf(name,"pixel_tot_avg_CS%i", i+1);
		saveToFile(config, this->name, "tot_analysis", name,
				h_pixel_tot_avg_CS[i]);
	}
	
	// alternative ToT pixel maps
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_sum_alt",
			h_pixel_tot_avg_alt);
	h_pixel_tot_avg_alt->Divide(h_pixel_tot_hits_alt);	// calculate average ToT
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_hits_alt",
			h_pixel_tot_hits_alt);
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_avg_alt",
			h_pixel_tot_avg_alt);

	// draw ToT pixel maps to file nicely
//	drawPixelNicely(config, "pixel_tot_avg_alt_nice", h_pixel_tot_avg_alt);

	// rebinned ToT pixel map	
	h_pixel_tot_avg_alt->Rebin2D(2,2);
	h_pixel_tot_avg_alt->Scale(1.0/(2*2));
	
//	drawPixelNicely(config, "pixel_tot_avg_alt_coarse_nice", h_pixel_tot_avg_alt);
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_avg_alt_coarse",
			h_pixel_tot_avg_alt);

	// 3D tot scatter histogram
	saveToFile(config, this->name, "tot_analysis", "pixel_tot_3d", h_pixel_tot_3d);
	
	// delete ToT maps
	delete h_module_tot_avg;
	delete h_module_tot_hits;
	delete h_pixel_tot_avg;
	delete h_pixel_tot_hits;
	delete h_pixel_tot_avg_alt;
	delete h_pixel_tot_hits_alt;
	
	for(int i=0; i<4; i++) {
		delete h_pixel_tot_hits_CS[i];
		delete h_pixel_tot_avg_CS[i];
	}
}

////////////////////////////////////////////////////////////////////////////////
//// CUT ANALYSIS //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// initialize
void Botho::init_cut_analysis(TbConfig &config) {
	/// control plots
	h_cuts_single = new TH1D("", ";;events", 7, 0, 7);
		h_cuts_single->GetXaxis()->SetBinLabel(1, "ALL EVENTS");
		h_cuts_single->GetXaxis()->SetBinLabel(2, "fTrack == kGood");
		h_cuts_single->GetXaxis()->SetBinLabel(3, "fClusters == kGood");
		h_cuts_single->GetXaxis()->SetBinLabel(4, "fTrackCentralRegion == kGood");
		h_cuts_single->GetXaxis()->SetBinLabel(5, "fTrackRegion == kGood");
		h_cuts_single->GetXaxis()->SetBinLabel(6, "fTrackMaskedRegion == kGood");
		h_cuts_single->GetXaxis()->SetBinLabel(7, "matched cluster");
	h_cuts_combined = new TH1D("", ";;events", 7, 0, 7);
		h_cuts_combined->GetXaxis()->SetBinLabel(1, "ALL EVENTS");
		h_cuts_combined->GetXaxis()->SetBinLabel(2, "fTrack == kGood");
		h_cuts_combined->GetXaxis()->SetBinLabel(3, "fClusters == kGood");
		h_cuts_combined->GetXaxis()->SetBinLabel(4, "fTrackCentralRegion == kGood");
		h_cuts_combined->GetXaxis()->SetBinLabel(5, "fTrackRegion == kGood");
		h_cuts_combined->GetXaxis()->SetBinLabel(6, "fTrackMaskedRegion == kGood");
		h_cuts_combined->GetXaxis()->SetBinLabel(7, "matched cluster");
		
	h_cuts_botho_single = new TH1D("", ";;events", 4, 0, 4);
		h_cuts_botho_single->GetXaxis()->SetBinLabel(1, "ALL EVENTS");
		h_cuts_botho_single->GetXaxis()->SetBinLabel(2, "fTrack == kGood");
		h_cuts_botho_single->GetXaxis()->SetBinLabel(3, "matched cluster");
		h_cuts_botho_single->GetXaxis()->SetBinLabel(4, "fTrackRegion == kGood");
	h_cuts_botho_combined = new TH1D("", ";;events", 4, 0, 4);
		h_cuts_botho_combined->GetXaxis()->SetBinLabel(1, "ALL EVENTS");
		h_cuts_botho_combined->GetXaxis()->SetBinLabel(2, "fTrack == kGood");
		h_cuts_botho_combined->GetXaxis()->SetBinLabel(3, "matched cluster");
		h_cuts_botho_combined->GetXaxis()->SetBinLabel(4, "fTrackRegion == kGood");
}		
	
/// analyse
void Botho::cut_analysis(const TbConfig &config, const Event &event) {

	// all events
	h_cuts_single->Fill("ALL EVENTS", 1);
	h_cuts_combined->Fill("ALL EVENTS", 1);
	h_cuts_botho_single->Fill("ALL EVENTS", 1);
	h_cuts_botho_combined->Fill("ALL EVENTS", 1);
	
	// single cuts
	if(event.fTrack == event::kGood) {
		h_cuts_single->Fill("fTrack == kGood", 1);
		h_cuts_botho_single->Fill("fTrack == kGood", 1);
	}
	if(event.fClusters == event::kGood) {
		h_cuts_single->Fill("fClusters == kGood", 1);
	}
	if(event.fTrackCentralRegion == event::kGood) {
		h_cuts_single->Fill("fTrackCentralRegion == kGood", 1);
	}
	if(event.fTrackRegion == event::kGood) {
		h_cuts_single->Fill("fTrackRegion == kGood", 1);
		h_cuts_botho_single->Fill("fTrackRegion == kGood", 1);
	}
	if(event.fTrackMaskedRegion == event::kGood) {
		h_cuts_single->Fill("fTrackMaskedRegion == kGood", 1);
	}
	if(cluster::getMatched(event.clusters, event) != -1) {
		h_cuts_single->Fill("matched cluster", 1);
		h_cuts_botho_single->Fill("matched cluster", 1);
	}

	// combined cuts
	if(event.fTrack == event::kGood) {
		h_cuts_combined->Fill("fTrack == kGood", 1);
		if(event.fClusters == event::kGood) {
			h_cuts_combined->Fill("fClusters == kGood", 1);
			if(event.fTrackCentralRegion == event::kGood) {
				h_cuts_combined->Fill("fTrackCentralRegion == kGood", 1);
				if(event.fTrackRegion == event::kGood) {
					h_cuts_combined->Fill("fTrackRegion == kGood", 1);
					if(event.fTrackMaskedRegion == event::kGood) {
						h_cuts_combined->Fill("fTrackMaskedRegion == kGood", 1);
						if(cluster::getMatched(event.clusters, event) != -1) {
							h_cuts_combined->Fill("matched cluster", 1);
						}
					}
				}
			}
		}
	}
	
	if(event.fTrack == event::kGood) {
		h_cuts_botho_combined->Fill("fTrack == kGood", 1);
		if(cluster::getMatched(event.clusters, event) != -1) {
			h_cuts_botho_combined->Fill("matched cluster", 1);
			if(event.fTrackRegion == event::kGood) {
				h_cuts_botho_combined->Fill("fTrackRegion == kGood", 1);
			}
		}
	}
}

/// finalize
void Botho::finalize_cut_analysis(const TbConfig &config) {
	// save histograms
	saveToFile(config, this->name, "cut_analysis", "cuts_single",
			h_cuts_single);
	saveToFile(config, this->name, "cut_analysis", "cuts_combined",
			h_cuts_combined);
	saveToFile(config, this->name, "cut_analysis", "cuts_botho_single",
			h_cuts_botho_single);
	saveToFile(config, this->name, "cut_analysis", "cuts_botho_combined",
			h_cuts_botho_combined);
	
	// free memory
	delete h_cuts_single;
	delete h_cuts_combined;
	delete h_cuts_botho_single;
	delete h_cuts_botho_combined;
}

////////////////////////////////////////////////////////////////////////////////
/// HIT EFFICIENCY ANALYSIS ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// initialize
void Botho::init_efficiency_analysis(TbConfig &config) {
	// get DUT dependend variables
	const DUT* dut = config.getDut(this->iden);
	// get pitches
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();
	// get column and row numbers
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();

	// reset efficiency variables
	n_tracks_raw = 0;
	n_tracks_good = 0;
	n_hits = 0;
	
	// get command line parameters
	do_eff_square_selection =
			(bool) config.cmdLineExtras_argGetter(
			(string) "A_botho_DoSquareSelection", (bool) false,
			(string) "switches on doing square area depending analysis");

	/// efficiency module map histograms
	h_trackModuleMap = new TH2D("", "Module Track Map; Column; Row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
	h_hitModuleMap = new TH2D("", "Module Hit Map; Column; Row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
			
	
	/// efficiency pixel map histograms
	h_trackPixelMap = new TH2D("",
			"Track Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
	h_hitPixelMap = new TH2D("",
			"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
			
	/// square selected pixel maps		
	if(do_eff_square_selection) {
		/// efficiency pixel map histograms
		h_trackPixelMap_square_selected = new TH2D("",
				"Track Pixel Map;Long side [#mum];Short Side[#mum]",
				TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
		h_hitPixelMap_square_selected = new TH2D("",
				"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
				TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
	}
			
	/// efficiency2.cc analysis rebuild ----------------------------------------
	eff2cc_ntracks = 0;
	eff2cc_nhits = 0;
	eff2cc_anyhit = 0;
	
	h_eff2cc_trackMap = new TH2D("", "Track Map;Column;Row", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_eff2cc_hitMap = new TH2D("", "Hit Map;Column;Row", n_cols, 0, n_cols, 
			n_rows, 0, n_rows);
	h_eff2cc_anyhitMap = new TH2D("", "Any Hit Map;Column;Row", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
			
	h_eff2cc_trackPixelMap = new TH2D("",
			"Track Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
	h_eff2cc_hitPixelMap = new TH2D("",
			"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
	
	if(do_eff_square_selection) {
		/// efficiency pixel map histograms
		h_eff2cc_trackPixelMap_square_selected = new TH2D("",
				"Track Pixel Map;Long side [#mum];Short Side[#mum]",
				TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
		h_eff2cc_hitPixelMap_square_selected = new TH2D("",
				"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
				TMath::Nint(pitchX), 0, pitchX, TMath::Nint(pitchY), 0, pitchY);
	}

}			

// analyze
void Botho::efficiency_analysis(const TbConfig &config, const Event &event) {
	
	////////////////////////////////////////////
	/// RAW DATA ANALYSIS                    ///
	////////////////////////////////////////////

/*	
//	------------------------------------------------------------
//  TEMPORARY STUFF to make botho efficiency analysis equal efficiency2.cc analysis
//	to quickly get the results for RD50 June 2014 talk on highly irrad. CIS modules
//Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		TBALOG ( kDEBUG3 ) << "Event not passed all cuts" << endl;
		return;
	}
	//Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		TBALOG ( kDEBUG3 ) << "Clusters not present in event" << endl;
		return;
	}
	// Are we in the central region on the chip?
	if (event.fTrackCentralRegion != event::kGood) {
		return;
	}
	// Are we in a good region of the chip?
	if (event.fTrackRegion != event::kGood) {
		return;
	}

	//avoid masked regions
	if (event.fTrackMaskedRegion != event::kGood) {
		return;
	}
//	------------------------------------------------------------
//  END OF TEMPORARY STUFF
*/	
	
	// increment number of raw tracks
	n_tracks_raw++;
	
	////////////////////////////////////////////
	/// CUT 1:  EVENT.FTRACK == EVENT::KGOOD ///
	////////////////////////////////////////////
	
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {  // if not event.fTrack == event::kGood
		return;                         // leave event
	}
	
	/////////////////////
	//// CUT 1 ANALYSIS -------------------------------------------------------
	/////////////////////
	
	// increase number of good tracks
	n_tracks_good++;
	
	// module track map
	h_trackModuleMap->Fill(tbutils::getCol(event.trackX, event),
			tbutils::getRow(event.trackY, event));
	
	// produce pixel track map
	h_trackPixelMap->Fill(tbutils::getFoldedX(event), tbutils::getPixelY(event));
	// square selected pixel track map
	if(do_eff_square_selection) {
		if(in_square(ROI_square,
				tbutils::getCol(event.trackX, event),
				tbutils::getRow(event.trackY, event))) {
			h_trackPixelMap_square_selected->
					Fill(tbutils::getFoldedX(event),
						tbutils::getPixelY(event));
		}
	}
	
	// efficiency2.cc rebuilt analysis
	//Check that clusters have been made successfully
	if((event.fClusters == event::kGood) &&
	// Are we in the central region on the chip?
	  (event.fTrackCentralRegion == event::kGood) &&
	// Are we in a good region of the chip?
//	  (event.fTrackRegion == event::kGood) &&   // replace with new geometric cut
	  ((this->*track_check_central_region)(event) == event::kGood) && 
	//avoid masked regions
	  (event.fTrackMaskedRegion == event::kGood))
	{
		// increase number of good tracks
		eff2cc_ntracks++;
	
		// module track map
		h_eff2cc_trackMap->Fill(tbutils::getCol(event.trackX, event),
				tbutils::getRow(event.trackY, event));
	
		// produce pixel track map
		h_eff2cc_trackPixelMap->Fill(tbutils::getFoldedX(event),
				tbutils::getPixelY(event));
		// square selected pixel track map
		if(do_eff_square_selection) {
			if(in_square(ROI_square,
					tbutils::getCol(event.trackX, event),
					tbutils::getRow(event.trackY, event))) {
				h_eff2cc_trackPixelMap_square_selected->
						Fill(tbutils::getFoldedX(event),
							tbutils::getPixelY(event));
			}
		}
		
		//Are there any hits in the DUT
		if (event.hits.size() > 0) {
			eff2cc_anyhit++;
			//TBALOG( kDEBUG2 ) << "Found hit in DUT!" << endl;
		}
		
		//plot the position of *any* single hits
		for (vector<PllHit*>::const_iterator i = event.hits.begin();
				i != event.hits.end(); ++i) {
			if ((*i)->iden != this->iden)
				continue;
			h_eff2cc_anyhitMap->Fill((*i)->col - 1, (*i)->row - 1); //minus 1 due to fill counting from 1
		}
	}
	
	////////////////////////////////////////////
	/// CUT 2:  MATCHED CLUSTER              ///
	////////////////////////////////////////////
	
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
	                       // of matched cluster in vector of clusters of event
	if (matched == -1) {                              
		return;	// leave event
	}
	
	/////////////////////
	//// CUT 2 ANALYSIS -------------------------------------------------------
	/////////////////////
	
	// increase number of matched hits
	n_hits++;
	
	// fill module hit map
	h_hitModuleMap->Fill(tbutils::getCol(event.trackX, event),
			tbutils::getRow(event.trackY, event));
	
	// fill pixel hit map
	h_hitPixelMap->Fill(tbutils::getFoldedX(event), tbutils::getPixelY(event));
	// square selected pixel hit map
	if(do_eff_square_selection) {
		if(in_square(ROI_square,
				tbutils::getCol(event.trackX, event),
				tbutils::getRow(event.trackY, event))) {
			h_hitPixelMap_square_selected->
					Fill(tbutils::getFoldedX(event),
						tbutils::getPixelY(event));
		}
	}
	
	// efficiency2.cc rebuilt analysis
	if((event.fClusters == event::kGood) &&
	  (event.fTrackCentralRegion == event::kGood) &&
	  (event.fTrackRegion == event::kGood) &&
	  (event.fTrackMaskedRegion == event::kGood))
	{
		// increase number of matched hits
		eff2cc_nhits++;
	
		// fill module hit map
		h_eff2cc_hitMap->Fill(tbutils::getCol(event.trackX, event),
				tbutils::getRow(event.trackY, event));
		// LEFT OUT -1 in hitMap that is used in original analysis
	
		// fill pixel hit map
		h_eff2cc_hitPixelMap->Fill(tbutils::getFoldedX(event),
				tbutils::getPixelY(event));
		// square selected pixel hit map
		if(do_eff_square_selection) {
			if(in_square(ROI_square,
					tbutils::getCol(event.trackX, event),
					tbutils::getRow(event.trackY, event))) {
				h_eff2cc_hitPixelMap_square_selected->
						Fill(tbutils::getFoldedX(event),
							tbutils::getPixelY(event));
			}
		}
	}
	
}

// finalize
void Botho::finalize_efficiency_analysis(const TbConfig &config) {
	// efficiencies
 	
 	TBALOG( kINFO ) << " raw racks:          " << n_tracks_raw << endl;
	TBALOG( kINFO ) << " good tracks:        " << n_tracks_good << endl;
	TBALOG( kINFO ) << " matched hits:       " << n_hits << endl;
	TBALOG( kINFO ) << " Efficiency (match): "
			<< double(n_hits) / double(n_tracks_good) << endl;
//	TBALOG( kINFO ) << " Efficiency (any):   "
//			<< double(anyhit) / double(ntracks) << endl;
 	
 	// efficiency module maps
	saveToFile(config, this->name, "efficiency", "trackModuleMap", h_trackModuleMap);
	saveToFile(config, this->name, "efficiency", "hitModuleMap", h_hitModuleMap);
	
	TH2D *h_effModuleMap = new TH2D("","Hit Efficiency Module Map;Column;Row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
	h_effModuleMap->Divide(h_hitModuleMap, h_trackModuleMap, 1.0, 1.0);
	saveToFile(config, this->name, "efficiency", "effModuleMap", h_effModuleMap);
	delete h_effModuleMap;
	
 	// delete efficiency module maps
	delete h_trackModuleMap;
	delete h_hitModuleMap;
 	
 	// rebin pixel maps
 	//double rebinX = 10;
 	double rebinX = 5; //(botho 2014-06-05)
 	double rebinY = 5;
 	
 	h_trackPixelMap->Rebin2D(rebinX,rebinY);
 	h_hitPixelMap->Rebin2D(rebinX,rebinY);
 	if(do_eff_square_selection) {
 		h_trackPixelMap_square_selected->Rebin2D(rebinX,rebinY);
	 	h_hitPixelMap_square_selected->Rebin2D(rebinX,rebinY);
	}
 	
 	
 	TH2D *h_effPixelMap = new TH2D("",
			"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(pitchX)/rebinX, 0, pitchX, TMath::Nint(pitchY)/rebinY, 0, pitchY);
	h_effPixelMap->Divide(h_hitPixelMap, h_trackPixelMap, 1.0, 1.0);
	
	// efficiency pixel maps
	saveToFile(config, this->name, "efficiency", "trackPixelMap", h_trackPixelMap);
	saveToFile(config, this->name, "efficiency", "hitPixelMap", h_hitPixelMap);
	saveToFile(config, this->name, "efficiency", "effPixelMap", h_effPixelMap);
	delete h_effPixelMap;
	
	// delete efficiency pixel maps
	delete h_trackPixelMap;
	delete h_hitPixelMap;
	
	if(do_eff_square_selection) {
		TH2D *h_effPixelMap_square_selected = new TH2D("",
				"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
				TMath::Nint(pitchX)/rebinX, 0,
				pitchX, TMath::Nint(pitchY)/rebinY, 0, pitchY);
		h_effPixelMap_square_selected->Divide(
				h_hitPixelMap_square_selected,
				h_trackPixelMap_square_selected,
				1.0, 1.0);
		// efficiency pixel maps
		saveToFile(config, this->name, "efficiency",
				"trackPixelMap_square_selected", 
				h_trackPixelMap_square_selected);
		saveToFile(config, this->name, "efficiency",
				"hitPixelMap_square_selected",
				h_hitPixelMap_square_selected);
		saveToFile(config, this->name, "efficiency",
				"effPixelMap_square_selected", 
				h_effPixelMap_square_selected);
		delete h_effPixelMap_square_selected;
	
		// delete efficiency pixel maps
		delete h_trackPixelMap_square_selected;
		delete h_hitPixelMap_square_selected;
	}
	
	
	///////////////////////
	//------------------------------------------------------------
	// efficienciy2.cc rebuild
 	
 	TBALOG( kINFO ) << " tracks:             " << eff2cc_ntracks << endl;
	TBALOG( kINFO ) << " tracks w/ hit:      " << eff2cc_nhits   << endl;
	TBALOG( kINFO ) << " hits:               " << eff2cc_anyhit  << endl;
	TBALOG( kINFO ) << " Efficiency (match): "
			<< double(eff2cc_nhits) / double(eff2cc_ntracks) << endl;
	TBALOG( kINFO ) << " Efficiency (any):   "
			<< double(eff2cc_anyhit) / double(eff2cc_ntracks) << endl;
 	
 	// efficiency module maps
	saveToFile(config, this->name, "efficiency", "eff2cc_trackModuleMap",
			h_eff2cc_trackMap);
	saveToFile(config, this->name, "efficiency", "eff2cc_hitModuleMap",
			h_eff2cc_hitMap);
	
	TH2D *h_eff2cc_effModuleMap = new TH2D("",
			"Hit Efficiency Module Map;Column;Row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
	h_eff2cc_effModuleMap->Divide(h_eff2cc_hitMap, h_eff2cc_trackMap, 1.0, 1.0);
	saveToFile(config, this->name, "efficiency", "eff2cc_effModuleMap",
			h_eff2cc_effModuleMap);
	delete h_eff2cc_effModuleMap;
	
 	// delete efficiency module maps
	delete h_eff2cc_trackMap;
	delete h_eff2cc_hitMap;
 	
 	// rebin pixel maps
 	//double rebinX = 10;
 	rebinX = 5; //(botho 2014-06-05)
 	rebinY = 5;
 	
 	h_eff2cc_trackPixelMap->Rebin2D(rebinX,rebinY);
 	h_eff2cc_hitPixelMap->Rebin2D(rebinX,rebinY);
 	if(do_eff_square_selection) {
 		h_eff2cc_trackPixelMap_square_selected->Rebin2D(rebinX,rebinY);
	 	h_eff2cc_hitPixelMap_square_selected->Rebin2D(rebinX,rebinY);
	}
 	
 	
 	TH2D *h_eff2cc_effPixelMap = new TH2D("",
			"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(pitchX)/rebinX, 0, pitchX, TMath::Nint(pitchY)/rebinY, 0, pitchY);
	h_eff2cc_effPixelMap->Divide(h_eff2cc_hitPixelMap,
			h_eff2cc_trackPixelMap, 1.0, 1.0);
	
	// efficiency pixel maps
	saveToFile(config, this->name, "efficiency", "eff2cc_trackPixelMap",
			h_eff2cc_trackPixelMap);
	saveToFile(config, this->name, "efficiency", "eff2cc_hitPixelMap",
			h_eff2cc_hitPixelMap);
	saveToFile(config, this->name, "efficiency", "eff2cc_effPixelMap",
			h_eff2cc_effPixelMap);
	delete h_eff2cc_effPixelMap;
	
	// delete efficiency pixel maps
	delete h_eff2cc_trackPixelMap;
	delete h_eff2cc_hitPixelMap;
	
	if(do_eff_square_selection) {
		TH2D *h_eff2cc_effPixelMap_square_selected = new TH2D("",
				"Hit Pixel Map;Long side [#mum];Short Side[#mum]",
				TMath::Nint(pitchX)/rebinX, 0,
				pitchX, TMath::Nint(pitchY)/rebinY, 0, pitchY);
		h_eff2cc_effPixelMap_square_selected->Divide(
				h_eff2cc_hitPixelMap_square_selected,
				h_eff2cc_trackPixelMap_square_selected,
				1.0, 1.0);
		// efficiency pixel maps
		saveToFile(config, this->name, "efficiency",
				"eff2cc_trackPixelMap_square_selected", 
				h_eff2cc_trackPixelMap_square_selected);
		saveToFile(config, this->name, "efficiency",
				"eff2cc_hitPixelMap_square_selected",
				h_eff2cc_hitPixelMap_square_selected);
		saveToFile(config, this->name, "efficiency",
				"eff2cc_effPixelMap_square_selected", 
				h_eff2cc_effPixelMap_square_selected);
		delete h_eff2cc_effPixelMap_square_selected;
	
		// delete efficiency pixel maps
		delete h_eff2cc_trackPixelMap_square_selected;
		delete h_eff2cc_hitPixelMap_square_selected;
	}
}	

////////////////////////////////////////////////////////////////////////////////
/// EDGE EFFICIENCY ANALYSIS ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// initialize edge efficiency analysis
void Botho::init_edgeefficiency_analysis(TbConfig &config) {
	const DUT *dut = config.getDut(this->iden);
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();
	epitchX = dut->getePitchX();
	epitchY = dut->getePitchY();
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();
	skip_cols = dut->getSkipCols();
	skip_rows = dut->getSkipRows();
	
	ee_ntracks = 0;
	ee_nTracksEdge = 0;
	ee_nHitsEdge = 0;
	
	// initialize histograms Stefano's analysis
	h_trackOutLeft = new TH1D("", "Hit Edge Pixel Map Left;Long side [#mum];",
			110, pitchX/2.0 - epitchX-500.0, pitchX/2.0);
	h_hitOutLeft = new TH1D("", "Hit Edge Pixel Map Left;Long side [#mum];",
			110, pitchX/2.0 - epitchX-500.0, pitchX/2.0);
	h_effOutLeft = new TH1D("", "Hit Edge Pixel Map Left;Long side [#mum];",
			110, pitchX/2.0 - epitchX-500.0, pitchX/2.0);
	h_trackOutRight = new TH1D("",
			"Hit Edge Pixel Map Right;Long side [#mum];", 110,
			(n_cols-1)*pitchX - pitchX/2.0,
			(n_cols-1)*pitchX+epitchX-pitchX/2.0 + 500.0);
	h_hitOutRight = new TH1D("",
			"Hit Edge Pixel Map Right;Long side [#mum];", 110,
			(n_cols-1)*pitchX - pitchX/2.0,
			(n_cols-1)*pitchX+epitchX-pitchX/2.0 + 500.0);
	h_effOutRight = new TH1D("",
			"Hit Edge Pixel Map Right;Long side [#mum];", 110,
			(n_cols-1)*pitchX - pitchX/2.0,
			(n_cols-1)*pitchX+epitchX-pitchX/2.0 + 500.0);
			
	// original edge analysis histograms
	h_EEAnyHitSensorMap = new TH2D("", "Any Edge Hit Map;Column;Row",
			n_cols, 0, n_cols, n_rows, 0, n_rows);
	h_EEHitSensorMap = new TH2D("", "Hit Map Edge Tracks;Column;Row",
			n_cols + 4, -2, n_cols + 2, n_rows, 0, n_rows);
	h_EETrackSensorMap = new TH2D("", "Track Map Edge Tracks;Column;Row",
			n_cols + 4, -2, n_cols + 2, n_rows, 0, n_rows);

/*	h_EESensorMap = new TH2D("", "Hit Map Edge Tracks;Column;Row",
			n_cols + 4, -2, n_cols + 2, n_rows, 0, n_rows);
*/

	h_EEHitPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EEHitLeftPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EEHitRightPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);

	h_EETrackPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EETrackLeftPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	h_EETrackRightPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), 0, epitchX + pitchX,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	// different binning for 125um module of September 2012 PPS TB at CERN
	if((strcmp(config.tbslot, "eudetPPSSep2012FEI3") == 0 && this->iden == 10)
		|| (strcmp(config.tbslot, "eudetPPSNov2013FEI34") == 0 && this->iden == 10)) {
		TBALOG (kINFO) << "VTT-E3 125 um edge in September 2012 run at CERN"
					<< "or November 2013 run at DESY"
					<< std::endl;
		delete h_EETrackRightPixelMap;
		delete h_EEHitRightPixelMap;
		h_EETrackRightPixelMap = new TH2D("",
			"Track Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), -2.5, epitchX + pitchX -2.5,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
		h_EEHitRightPixelMap = new TH2D("",
			"Hit Edge Pixel Map;Long side [#mum];Short Side[#mum]",
			TMath::Nint(epitchX + pitchX), -2.5, epitchX + pitchX -2.5,
			TMath::Nint(pitchY) * 2, 0, pitchY * 2);
	}
			
	// initialize histograms Botho's analysis
	// number of bins
	const int bothoEE_nbins = 220; // gives bin size of 5 um
	h_bothoEE_trackOutLeft = new TH1D("", "Hit Edge Pixel Map Left;Long side [#mum];",
			bothoEE_nbins, - epitchX-100.0, pitchX);
	h_bothoEE_hitOutLeft = new TH1D("", "Hit Edge Pixel Map Left;Long side [#mum];",
			bothoEE_nbins, - epitchX-100.0, pitchX);
	h_bothoEE_effOutLeft = new TH1D("", "Hit Edge Pixel Map Left;Long side [#mum];",
			bothoEE_nbins, - epitchX-100.0, pitchX);
	h_bothoEE_trackOutRight = new TH1D("",
			"Hit Edge Pixel Map Right;Long side [#mum];", bothoEE_nbins,
			(n_cols-2)*pitchX,
			(n_cols-1)*pitchX+epitchX +100.0);
	h_bothoEE_hitOutRight = new TH1D("",
			"Hit Edge Pixel Map Right;Long side [#mum];", bothoEE_nbins,
			(n_cols-2)*pitchX,
			(n_cols-1)*pitchX+epitchX +100.0);
	h_bothoEE_effOutRight = new TH1D("",
			"Hit Edge Pixel Map Right;Long side [#mum];", bothoEE_nbins,
			(n_cols-2)*pitchX,
			(n_cols-1)*pitchX+epitchX +100.0);
}

// do edge efficiency analysis for event
void Botho::edgeefficiency_analysis(const TbConfig &config, const Event &event) {
	
	// check that track has passed all cuts
	if(event.fTrack != event::kGood)
		return;
		
	ee_ntracks++;
	
/*	// Are we *not* in the central region on the chip? 
	if (event.fTrackCentralRegion == event::kGood) {
		TBALOG (kDEBUG3) << "Track in central region of chip" << endl;
		return;
	}
*/

	// The edge in Y (ganged- and inter-ganged pixels) should still be ignored
/*	if ((event.trackY + 0.5 * pitchY) < (pitchY * skip_rows)) {
		return;
	}
	if ((event.trackY + 0.5 * pitchY) > (pitchY * (n_rows - skip_rows))) {
		return;
	}
*/
	// new boundaries in Y, just 7 rows on top because they are ganged,
	// and the bottom row, because there are no hits from the reconstruction
	// (changed 2014-07-24)
	if ((event.trackY + 0.5 * pitchY) < (pitchY * 1)) {
		return;
	}
	if ((event.trackY + 0.5 * pitchY) > (pitchY * (n_rows-7)/*(n_rows - 8)*/)) {
		return;
	}
/*	if(event.dut->getDUTid() == 10) {
		if ((event.trackY + 0.5 * pitchY) < (pitchY * 30)) {
			return;
		}
		if ((event.trackY + 0.5 * pitchY) > (pitchY * (153))) {
			return;
		}
	}
*/	
	// Stefano's version -------------------------------------------------
	float col;
	if( ( pitchX/2.-epitchX-500 < event.trackX ) && ( event.trackX < pitchX/2. ) )
	{
		h_trackOutLeft->Fill(event.trackX);
		col = 0;
		bool found = false;
		for( int ii = 0; ii < event.clusters.size(); ii++)
		{
			for(int cc = 0; cc < event.clusters.at( ii ).size(); cc++)
			{
				if( event.clusters.at( ii ).at( cc )->col == col )
				{
					if( cluster::getMatched( event.clusters, event ) != -1 )
					{
					   	h_hitOutLeft->Fill(event.trackX);
						found = true;
					   	break;
					}
				}
			}
			if(found) break;
		}
	}
	else
	// select tracks for RIGHT edge of the sensor
	if( ( (n_cols-1)*pitchX-pitchX/2. < event.trackX ) &&
			( event.trackX < (n_cols-1)*pitchX+epitchX-pitchX/2.+500. ) )
	{
		h_trackOutRight->Fill(event.trackX);
		col = (n_cols-1);
		//col = 79;
		bool found = false;
		// Look for matched clusters only
		for( int ii = 0; ii < event.clusters.size(); ii++)
		{
			for(int cc = 0; cc < event.clusters.at( ii ).size(); cc++)
			{
				if( event.clusters.at( ii ).at( cc )->col == col ) //|| event.clusters.at( ii ).at( cc )->col == (col+1) )
				{
					if( cluster::getMatched( event.clusters, event ) != -1 )
					{
						h_hitOutRight->Fill(event.trackX);
						found = true;
						break;
					}
				}
			}
			if(found) break;
		}
	}  
	
	// Botho's version -------------------------------------------------------
	//float col;
	if( (-epitchX-100 < event.trackX ) && ( event.trackX < pitchX ) )
	{
		h_bothoEE_trackOutLeft->Fill(event.trackX);
		col = 0;
		bool found = false;
		for( int ii = 0; ii < event.clusters.size(); ii++)
		{
			for(int cc = 0; cc < event.clusters.at( ii ).size(); cc++)
			{
				if( (event.clusters.at( ii ).at( cc )->col == col)
					||  (event.clusters.at( ii ).at( cc )->col == col+1))
				{
					if( cluster::getMatched( event.clusters, event ) != -1 )
					{
					   	h_bothoEE_hitOutLeft->Fill(event.trackX);
						found = true;
					   	break;
					}
				}
			}
			if(found) break;
		}
	}
	else
	// select tracks for RIGHT edge of the sensor
	if( ( (n_cols-2)*pitchX < event.trackX ) &&
			( event.trackX < (n_cols-1)*pitchX+epitchX+100. ) )
	{
		h_bothoEE_trackOutRight->Fill(event.trackX);
		col = (n_cols-1);
		//col = 79;
		bool found = false;
		// Look for matched clusters only
		for( int ii = 0; ii < event.clusters.size(); ii++)
		{
			for(int cc = 0; cc < event.clusters.at( ii ).size(); cc++)
			{
				if( (event.clusters.at( ii ).at( cc )->col == col )
					|| ( event.clusters.at( ii ).at( cc )->col == col-1) ) //|| event.clusters.at( ii ).at( cc )->col == (col+1) )
				{
					if( cluster::getMatched( event.clusters, event ) != -1 )
					{
						h_bothoEE_hitOutRight->Fill(event.trackX);
						found = true;
						break;
					}
				}
			}
			if(found) break;
		}
	}
	
//	--------------------------------------------------------------------------
//	ORIGINAL VERSION FROM TRUNK EDGEEFFICIENCY.CC (2014-06-17)
	
	// Make sure we are not too far off the sensor in the long (X) direction for edge studies
	// (for the X edge we want to include half a pixel out)  
	// select tracks for LEFT edge of the sensor
	if ((-1.0 * epitchX < event.trackX) && (event.trackX < pitchX)) {
		h_EETrackPixelMap->Fill(event.trackX + epitchX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));
		h_EETrackLeftPixelMap->Fill(event.trackX + epitchX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));

		// Look for matched clusters only
		if (cluster::getMatched(event.clusters, event) != -1) {
			TBALOG(kDEBUG2) << "Found matching cluster!" << endl;
			h_EEHitPixelMap->Fill(event.trackX + epitchX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
			h_EEHitLeftPixelMap->Fill(event.trackX + epitchX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
		}
	} else
	// select tracks for RIGHT edge of the sensor
	if (((n_cols - 2) * pitchX < event.trackX)
			&& (event.trackX < (n_cols - 1) * pitchX + epitchX)) {
		h_EETrackPixelMap->Fill(((n_cols - 1) * pitchX + epitchX) - event.trackX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));
		h_EETrackRightPixelMap->Fill(event.trackX - (n_cols - 2) * pitchX,
				tbutils::getTrackTwoPixelEdgeDistance(event.trackY + epitchY,
						pitchY));

		// Look for matched clusters only
		if (cluster::getMatched(event.clusters, event) != -1) {
			TBALOG(kDEBUG2) << "Found matching cluster!" << endl;
			h_EEHitPixelMap->Fill(
					((n_cols - 1) * pitchX + epitchX) - event.trackX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
			h_EEHitRightPixelMap->Fill(event.trackX - (n_cols - 2) * pitchX,
					tbutils::getTrackTwoPixelEdgeDistance(
							event.trackY + epitchY, pitchY));
		}
	} else {
		return;
	}

	// Fill Sensor Maps using both sides of the sensor
	// Count the tracks
	ee_nTracksEdge++;

	// plot the position of *any* single hits
	for (vector<PllHit*>::const_iterator i = event.hits.begin();
			i != event.hits.end(); ++i) {
		if ((*i)->iden != this->iden)
			continue;
		h_EEAnyHitSensorMap->Fill((*i)->col, (*i)->row);
	}

	col = (event.trackX + 0.5 * pitchX) / pitchX;
	float row = (event.trackY + 0.5 * pitchY) / pitchY;

	// Fill the Edge Efficiency Sensor Map
	h_EETrackSensorMap->Fill(col, row);

	// Are there any hits in the DUT?
	if (event.hits.size() > 0) {
		TBALOG(kDEBUG2) << "Found hit in DUT!" << endl;
	}

	// Look for matched clusters
	if (cluster::getMatched(event.clusters, event) == -1) {
		return;
	}
	TBALOG(kDEBUG2) << "Found matching cluster!" << endl;

	// count the hits 
	ee_nHitsEdge++;
	h_EEHitSensorMap->Fill(col, row);
//	h_EESensorMap->SetStats(false);
//	h_EESensorMap->Divide(h_EEHitSensorMap, h_EETrackSensorMap, 1.0, 1.0, "B");
}

// finalize edge efficiency analysis
void Botho::finalize_edgeefficiency_analysis(const TbConfig &config) {

	TBALOG(kINFO) << " Tracks:             " << ee_ntracks << endl;
	TBALOG(kINFO) << " Edge tracks:        " << ee_nTracksEdge << endl;
	TBALOG(kINFO) << " Edge tracks w/ hit: " << ee_nHitsEdge << endl;
	
	// calculate Stefano's edge efficiencies ----------------------------------
	h_effOutRight->Divide(h_hitOutRight, h_trackOutRight, 1.0,1.0, "B");
	h_effOutLeft->Divide(h_hitOutLeft, h_trackOutLeft, 1.0, 1.0, "B");	
	
	// save histograms
	saveToFile(config, this->name, "edgeefficiency", "trackOutLeft_fine",
			h_trackOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "hitOutLeft_fine",
			h_hitOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "effOutLeft_fine",
			h_effOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "trackOutRight_fine",
			h_trackOutRight);
	saveToFile(config, this->name, "edgeefficiency", "hitOutRight_fine",
			h_hitOutRight);
	saveToFile(config, this->name, "edgeefficiency", "effOutRight_fine",
			h_effOutRight);
	
	// REBIN	
	h_trackOutLeft->Rebin(5);
	h_trackOutLeft->Scale(1.0/5.0);
	h_hitOutLeft->Rebin(5);
	h_hitOutLeft->Scale(1.0/5.0);
	h_effOutLeft->Rebin(5);
	h_effOutLeft->Scale(1.0/5.0);
	
	h_trackOutRight->Rebin(5);
	h_trackOutRight->Scale(1.0/5.0);
	h_hitOutRight->Rebin(5);
	h_hitOutRight->Scale(1.0/5.0);
	h_effOutRight->Rebin(5);
	h_effOutRight->Scale(1.0/5.0);
	
	// save coarse histograms
	saveToFile(config, this->name, "edgeefficiency", "trackOutLeft",
			h_trackOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "hitOutLeft",
			h_hitOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "effOutLeft",
			h_effOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "trackOutRight",
			h_trackOutRight);
	saveToFile(config, this->name, "edgeefficiency", "hitOutRight",
			h_hitOutRight);
	saveToFile(config, this->name, "edgeefficiency", "effOutRight",
			h_effOutRight);
			
	// delete histograms
	delete h_trackOutLeft;
	delete h_hitOutLeft;
	delete h_effOutLeft;
	delete h_trackOutRight;
	delete h_hitOutRight;
	delete h_effOutRight;
	
	// calculate Botho's edge efficiencies ----------------------------------
	h_bothoEE_effOutRight->Divide(h_bothoEE_hitOutRight, h_bothoEE_trackOutRight, 1.0,1.0, "B");
	h_bothoEE_effOutLeft->Divide(h_bothoEE_hitOutLeft, h_bothoEE_trackOutLeft, 1.0, 1.0, "B");	
	
	// save histograms
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_trackOutLeft_fine",
			h_bothoEE_trackOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_hitOutLeft_fine",
			h_bothoEE_hitOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_effOutLeft_fine",
			h_bothoEE_effOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_trackOutRight_fine",
			h_bothoEE_trackOutRight);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_hitOutRight_fine",
			h_bothoEE_hitOutRight);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_effOutRight_fine",
			h_bothoEE_effOutRight);
	
	// REBIN	
	h_bothoEE_trackOutLeft->Rebin(5);
	h_bothoEE_trackOutLeft->Scale(1.0/5.0);
	h_bothoEE_hitOutLeft->Rebin(5);
	h_bothoEE_hitOutLeft->Scale(1.0/5.0);
	h_bothoEE_effOutLeft->Rebin(5);
	h_bothoEE_effOutLeft->Scale(1.0/5.0);
	
	h_bothoEE_trackOutRight->Rebin(5);
	h_bothoEE_trackOutRight->Scale(1.0/5.0);
	h_bothoEE_hitOutRight->Rebin(5);
	h_bothoEE_hitOutRight->Scale(1.0/5.0);
	h_bothoEE_effOutRight->Rebin(5);
	h_bothoEE_effOutRight->Scale(1.0/5.0);
	
	// save coarse histograms
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_trackOutLeft",
			h_bothoEE_trackOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_hitOutLeft",
			h_bothoEE_hitOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_effOutLeft",
			h_bothoEE_effOutLeft);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_trackOutRight",
			h_bothoEE_trackOutRight);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_hitOutRight",
			h_bothoEE_hitOutRight);
	saveToFile(config, this->name, "edgeefficiency", "bothoEE_effOutRight",
			h_bothoEE_effOutRight);
			
	// delete histograms
	delete h_bothoEE_trackOutLeft;
	delete h_bothoEE_hitOutLeft;
	delete h_bothoEE_effOutLeft;
	delete h_bothoEE_trackOutRight;
	delete h_bothoEE_hitOutRight;
	delete h_bothoEE_effOutRight;
	
	// original edgeefficiency.cc finalization --------------------------------
	
	// sensor maps
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeAnyHitSensorMap",
			h_EEAnyHitSensorMap);
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeTrackSensorMap",
			h_EETrackSensorMap);
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeHitSensorMap",
			h_EEHitSensorMap);
	// sensor efficiency
	TH2D *h_EESensorMap = (TH2D*) h_EEHitSensorMap->Clone("");
	h_EESensorMap->Divide(h_EEHitSensorMap, h_EETrackSensorMap, 1.0, 1.0, "B");
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffSensorMap",
			h_EESensorMap);
			
	// pixel maps
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeTrackPixelMap",
			h_EETrackPixelMap);
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeHitPixelMap",
			h_EEHitPixelMap);
	// edge efficiency
	h_EETrackPixelMap->Rebin2D(5,2);
	h_EEHitPixelMap->Rebin2D(5,2);
	TH2D *h_EEPixelMap = (TH2D*) h_EEHitPixelMap->Clone("");
	h_EEPixelMap->Divide(h_EEHitPixelMap, h_EETrackPixelMap, 1., 1., "B");
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffPixelMap",
			h_EEPixelMap);
	
	// left edge pixel maps
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeTrackLeftPixelMap",
			h_EETrackLeftPixelMap);
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeHitLeftPixelMap",
			h_EEHitLeftPixelMap);
	// edge efficiency
	h_EETrackLeftPixelMap->Rebin2D(5,2);
	h_EEHitLeftPixelMap->Rebin2D(5,2);
	TH2D *h_EELeftPixelMap = (TH2D*) h_EEHitLeftPixelMap->Clone("");
	h_EELeftPixelMap->Divide(h_EEHitLeftPixelMap, h_EETrackLeftPixelMap, 1., 1., "B");
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffLeftPixelMap",
			h_EELeftPixelMap);
			
	// right edge pixel maps
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeTrackRightPixelMap",
			h_EETrackRightPixelMap);
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeHitRightPixelMap",
			h_EEHitRightPixelMap);
	
	// edge efficiency
	h_EETrackRightPixelMap->Rebin2D(5,2);
	h_EEHitRightPixelMap->Rebin2D(5,2);
	
	TH2D *h_EERightPixelMap = (TH2D*) h_EEHitRightPixelMap->Clone("");
	h_EERightPixelMap->Divide(h_EEHitRightPixelMap, h_EETrackRightPixelMap, 1., 1., "B");
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffRightPixelMap",
			h_EERightPixelMap);
			
	/// PROJECTIONS		
	// both sides projection
	TH1D *h_EETrackPixelMapX = (TH1D*) h_EETrackPixelMap->ProjectionX(
			"TrackPixelMapX", -1, 999999, "E");
	TH1D *h_EEHitPixelMapX = (TH1D*) h_EEHitPixelMap->ProjectionX(
			"HitPixelMapX", -1, 999999, "E");
	TGraphAsymmErrors* tg_EEPixelMapX = new TGraphAsymmErrors(
			h_EEHitPixelMapX, h_EETrackPixelMapX,
			"cl=0.683 b(1,1) mode");
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffProjX",
			tg_EEPixelMapX);
					
	// left side projection
	TH1D *h_EETrackLeftPixelMapX = (TH1D*) h_EETrackLeftPixelMap->ProjectionX(
			"TrackLeftPixelMapX", -1, 999999, "E");
	TH1D *h_EEHitLeftPixelMapX = (TH1D*) h_EEHitLeftPixelMap->ProjectionX(
			"HitLeftPixelMapX", -1, 999999, "E");
	TGraphAsymmErrors* tg_EELeftPixelMapX = new TGraphAsymmErrors(
			h_EEHitLeftPixelMapX, h_EETrackLeftPixelMapX,
			"cl=0.683 b(1,1) mode");
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffLeftProjX",
			tg_EELeftPixelMapX);
			
	// right side projection
	TH1D *h_EETrackRightPixelMapX = (TH1D*) h_EETrackRightPixelMap->ProjectionX(
			"TrackRightPixelMapX", -1, 999999, "E");
	TH1D *h_EEHitRightPixelMapX = (TH1D*) h_EEHitRightPixelMap->ProjectionX(
			"HitRightPixelMapX", -1, 999999, "E");
	TGraphAsymmErrors* tg_EERightPixelMapX = new TGraphAsymmErrors(
			h_EEHitRightPixelMapX, h_EETrackRightPixelMapX,
			"cl=0.683 b(1,1) mode");
/*	// delete x errors of graph
	for(int i=0; i < tg_EERightPixelMapX->GetN(); i++) {
		tg_EERightPixelMapX->SetPointEXlow(i, 0);
		tg_EERightPixelMapX->SetPointEXhigh(i, 0);
	}
*/
	saveToFile(config, this->name, "edgeefficiency", "EE_edgeEffRightProjX",
			tg_EERightPixelMapX);
	
	// calculate asymmetric errors with calc_efficiency function
	TGraphAsymmErrors *tg_calceff_right = new TGraphAsymmErrors();
//				tg_EERightPixelMapX->Clone("calceff_right");
	for(int i = 1; i <= h_EETrackRightPixelMapX->GetNbinsX(); i++) {
		int ntracks_bin = h_EETrackRightPixelMapX->GetBinContent(i);
		int nhits_bin = h_EEHitRightPixelMapX->GetBinContent(i);
		tg_calceff_right->SetPoint(i,
				h_EETrackRightPixelMapX->GetBinCenter(i),
				(double)nhits_bin/ntracks_bin);
		double err_up, err_down;
		calc_efficiency(ntracks_bin, nhits_bin, err_up, err_down);
		tg_calceff_right->SetPointError(i, 2.5, 2.5, err_down, err_up);
	}
	saveToFile(config, this->name, "edgeefficiency", "EE_calceff_right",
			tg_calceff_right);
	delete tg_calceff_right;
	
	// delete histos
//	delete fileName;
	delete h_EETrackSensorMap;
	delete h_EEAnyHitSensorMap;
	delete h_EEHitSensorMap;
	delete h_EESensorMap;

	delete h_EETrackPixelMap;
	delete h_EETrackLeftPixelMap;
	delete h_EETrackRightPixelMap;

	delete h_EEHitPixelMap;
	delete h_EEHitLeftPixelMap;
	delete h_EEHitRightPixelMap;

	delete h_EEPixelMap;
	delete h_EELeftPixelMap;
	delete h_EERightPixelMap;

	delete h_EEHitPixelMapX;
	delete h_EEHitLeftPixelMapX;
	delete h_EEHitRightPixelMapX;

	delete h_EETrackPixelMapX;
	delete h_EETrackLeftPixelMapX;
	delete h_EETrackRightPixelMapX;

	delete tg_EEPixelMapX;
	delete tg_EELeftPixelMapX;
	delete tg_EERightPixelMapX;
}

// additional functions for edge analysis:

// error calculation copied from program by Lee Holloway 1989 provided by Stefano
void Botho::calc_efficiency (double Ntot, double Ngood, double &sup, double &sd) {

  // internal variables
  double term = 0.;
  double p = 0.;
  double d = 0.;
  double r = 0.;
  double s = 0.;
  double ans[2], tmp= 0.;
  int lim[2], i, j, k;

  // Set the value of b (choose among the following cases
  // if you want to change the confidence level used to
  // quote results):

   float b[2] = {.8415, .1585};   // for 68% C.L.
  //  float b[2] = {.975 , .025 };   // for 95% C.L.


  // Set up limits and check for special cases
  lim[0] = Ngood-1;
  lim[1] = Ngood;

 
  for (k = 0; k < 2; k++) { 
    if (Ngood == 0 && k==0) {
      ans[0] = 0;
      continue;
    }
    if (Ngood == Ntot && k==1) {
      ans[1] = 1;
      continue;
    }

    // 
    // Start a binary search procedure
    d=1.;
    p=1.;
    s=0.;    
    for (j = 0; j < 12; j++){
      cout << j << endl;
      d = TMath::Sign((d/2), (s-b[k]));
      p = p+d;      
      r = TMath::Log(p /(1.-p));
      term = Ntot * TMath::Log(1.-p);
      s = TMath::Exp(term);
      
      for (i = 1; i <= lim[k]; i++){
// 	cout << i << endl;
	tmp = (Ntot+1. - i)/i;
	term = term + r + TMath::Log(tmp);
	s = s + TMath::Exp(term);
      } // end loop on i

    } // end loop on j
    ans[k] = p;
        
  } // end loop on k
  
  double lower = ans[0];
  double upper = ans[1];
  

   //The following is done in the main routine calc_eff().
  
   double mean= Ngood / Ntot;  
   sup=upper-mean;
   sd=mean-lower;
/*   cout << Ngood << " " << Ntot << " " << upper << " " << lower << endl;
   cout << " mean = " << mean << " + " << sup << " - " << sd << endl;
*/

}


// copied from trunk version of edgeefficiency.cc (and slightly altered)
void Botho::drawToFileEdgeEff(const TbConfig& config,
		const char* analysisName, const char* histoName, TObject* histo,
		double fitLeft, double fitRight) const {

	const DUT* dut = config.getDut(this->iden);
	char* fileName = config.buildHistName(analysisName, histoName);
	double pitchX = dut->getPitchX();
	double epitchX = dut->getePitchX();

	// Make sure both a graph and hist can be used in the fit
	TGraphAsymmErrors* gr_histo = 0;
	TH1D* h_histo = 0;
	bool isGraph = false;
	if (histo->InheritsFrom("TGraph")) {
		isGraph = true;
		gr_histo = (TGraphAsymmErrors*) histo;
	} else if (histo->InheritsFrom("TH1D")) {
		h_histo = (TH1D*) histo;
	} else {
		TBALOG(kERROR)
				<< "You must have either a TH1D or TGraph for the edge fit!"
				<< endl;
		exit(1);
	}
	assert(tbutils::XOR(gr_histo, h_histo));

	TCanvas* canvas = new TCanvas("", "", 1200, 300);
	canvas->cd();
	gStyle->SetPadBottomMargin(0.15);
	gROOT->ForceStyle();
	TString title = ";Long pixel [#mum ]; Edge efficiency";
	TH1* h_plot = 0;
	if (isGraph) {
		gr_histo->SetLineWidth((Width_t) 1.3); // BJD: Really?!
		gr_histo->SetMarkerStyle(20);
		gr_histo->SetMarkerSize(0.7);
		gr_histo->Draw("AP");
		gr_histo->SetTitle(title);
		gr_histo->GetXaxis()->SetLabelSize(0.07);
		gr_histo->GetYaxis()->SetLabelSize(0.07);
		gr_histo->GetYaxis()->SetTitleSize(0.08);
		gr_histo->GetXaxis()->SetTitleSize(0.08);
		gr_histo->GetYaxis()->SetTitleOffset(0.4);
		gr_histo->GetXaxis()->SetTitleOffset(0.9);
		gr_histo->GetXaxis()->SetRangeUser(0., pitchX + epitchX);
		h_plot = gr_histo->GetHistogram();
	} else {
		h_histo->SetLineWidth((Width_t) 1.3); // BJD: Really?!
		h_histo->SetTitle(title);
		h_histo->GetYaxis()->SetTitleSize(0.09);
		h_histo->GetYaxis()->SetTitleOffset(0.04);
		h_histo->Draw();
		h_plot = h_histo;
	}

	// Draw an area covering the next to last pixel for visualization
	TH1D* h_innerPixel = (TH1D*) h_plot->Clone("innerPixel");
	for (int iBin = 1; iBin != h_innerPixel->GetNbinsX() + 1; ++iBin) {
		if ((iBin <= h_plot->FindBin(pitchX / 2))
				|| (iBin >= h_plot->FindBin(pitchX / 2 + epitchX)))
			h_innerPixel->SetBinContent(iBin, 0.0);
		else
			h_innerPixel->SetBinContent(iBin, h_innerPixel->GetMaximum());
	}
	//h_innerPixel->SetFillColor(kBlue);
	//h_innerPixel->SetLineColor(kGray);
	//h_innerPixel->SetFillStyle(3004);
	h_innerPixel->SetLineWidth((Width_t) 0.5); // BJD: Really?!
	h_innerPixel->Draw("same");

	TF1* fit = tbutils::GetGausLineFit(histo, "turnon_fullrange", fitLeft,
			fitRight, 50.0, 9.0, 0.98);
	fit->SetLineColor(kRed);
	fit->SetLineWidth((Width_t) 2.3); // BJD: Really?!
	fit->Draw("L same");
	double xf90 = fit->GetX(0.9);
	double xf10 = fit->GetX(0.1);
	TPad *pad2 = new TPad("", "", 0.405, 0.27, 0.705, 0.52);
	pad2->Draw();
	pad2->cd();

	TPaveStats *paramBox = new TPaveStats(0.0, 0.0, 1.0, 1.0);
	paramBox->AddText("Active edge fit");
	char * line = new char[300];
	sprintf(line, "Pixel width=%.1f#pm%.1f #mum", pitchX - fit->GetParameter(0),
			fit->GetParError(0));
	paramBox->AddText(line);
	sprintf(line, "Resolution=%.1f#pm%.1f #mum", fit->GetParameter(1),
			fit->GetParError(1));
	paramBox->AddText(line);
	sprintf(line, "Plateau eff.=%.3f#pm%.3f", fit->GetParameter(2),
			fit->GetParError(2));
	paramBox->AddText(line);
	sprintf(line, "90%% (10%%)= %.1f (%.1f) #mum", pitchX - xf90,
			pitchX - xf10);
	paramBox->AddText(line);
	paramBox->SetFillColor(0);
	paramBox->SetTextColor(kBlack);
	paramBox->Draw();
//	canvas->SaveAs(fileName);
	
	// added (2014-06-17) by Botho
	//config.saveToFile(this->name, histoName, (TNamed*) histo);
	config.saveToFile(this->name, histoName, canvas);

	delete canvas;
	delete fileName;
	delete fit;
	delete paramBox;
	return;

}

///////////////////////////////////////
/// PRIVATE MANAGEMENT FUNCTIONS ------
///////////////////////////////////////

void Botho::saveToFile(const TbConfig &config, const char *analysisName,
			const char *subfolderName, const char *histoName, TNamed *histo) {
	// basic function copied from config::saveToFile		
	config.tfile->cd();
	char* fullName = new char[500];
	if (config.organizeOutput) {
		config.tfile->mkdir(analysisName);
		config.tfile->cd(analysisName);
		sprintf(fullName, "%s", histoName);
	} else {
		sprintf(fullName, "%s-%s", analysisName, histoName);
	}
	for (int ii = 0; ii < 500; ii++) { // Histo not callable from the cint cli if '-' in the name
		if (fullName[ii] == '-') {
			fullName[ii] = '_';
		}
		if (fullName[ii] == '\0') {
			break;
		}
	}
	// new code by Botho
	// this creates the subfolder, if not yet existing, and makes it the
	// current working directory for the histogram to be saved
	gDirectory->mkdir(subfolderName);
	gDirectory->cd(subfolderName);
	// end of new code by Botho
	
	histo->SetName(fullName);
	histo->Write(fullName);
	delete[] fullName;
}

bool Botho::in_square(int *square, int x, int y) {
	// square must be integer array of format
	// int square[4] = {x1, y1, x2, y2};
	if (((x >= square[0]) && (x <= square[2]))
			&& ((y >= square[1]) && (y <= square[3]))) {
		return true;
	} else {
		return false;
	}
}

void Botho::drawPixelNicely(const TbConfig &config, const char *histoname, TH1 *histo) {
	// draw ToT pixel maps to file nicely
	// gStyle->Reset();
	Int_t gstyle_opt_stat = gStyle->GetOptStat();
   	
   	Bool_t groot_batchmode = gROOT->IsBatch();
   	gROOT->SetBatch();
   	gStyle->SetOptStat(0);
   
   	// create canvas
   	Double_t canvas_width = 900; // 1200
	Double_t canvas_height = 150;	// 200
	TCanvas *c2 = new TCanvas("qeff", "qeff", canvas_width, canvas_height);
	Double_t offset_width = canvas_width - c2->GetWw();
	Double_t offset_height = canvas_height - c2->GetWh();
	c2->SetWindowSize(canvas_width + offset_width, canvas_height + offset_height);	// make sure, actual canvas size is 800x600 (not important)
	c2->SetHighLightColor(2);
	c2->Range(0,0,1,1);
	c2->SetFillColor(0);
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
   
	// set margins
	c2->SetLeftMargin(0.05);
	c2->SetRightMargin(0.08);
	c2->SetTopMargin(0.05);
	c2->SetBottomMargin(0.25);
   
	c2->SetTickx(0);
	c2->SetTicky(0);

	//c2->SetFixedAspectRatio();
	Double_t asp_ratio = 6; // c2->GetAspectRatio();
   
   	// set line widths
	gStyle->SetLineWidth(2);
	c2->SetFrameLineWidth(2);

//	histo->GetXaxis()->SetTitle("Track x [#mum]");
	histo->GetXaxis()->SetNdivisions(405);
	histo->GetXaxis()->SetTickLength(-0.04);
	histo->GetXaxis()->SetLabelOffset(0.03);
	histo->GetXaxis()->SetLabelSize(0.12);
	histo->GetXaxis()->SetTitleSize(0.12);
	histo->GetXaxis()->SetTitleOffset(1.0);
//	histo->GetYaxis()->SetTitle("Track y [#mum]");
	histo->GetYaxis()->SetNdivisions(2,0,0,kFALSE);	// no "optimization"
	histo->GetYaxis()->SetLabelOffset(0.03/asp_ratio);
	histo->GetYaxis()->SetLabelSize(0.12);
	histo->GetYaxis()->SetTitleSize(0.12);
	histo->GetYaxis()->SetTitleOffset(0.20);
	histo->GetYaxis()->SetTickLength(histo->GetXaxis()->GetTickLength()*(1.0/asp_ratio));
//	histo->GetZaxis()->SetTitle("ToT");
	histo->GetZaxis()->SetNdivisions(5);

	histo->GetZaxis()->SetLabelOffset(0.006);
	histo->GetZaxis()->SetLabelSize(0.10);
	histo->GetZaxis()->SetTitleSize(0.12);
	histo->GetZaxis()->SetTitleOffset(0.22);
	histo->GetZaxis()->SetTickLength(-0.005);
	//   histo->GetZaxis()->SetRangeUser(8,26);  
	histo->Draw("COLZ");

	gPad->Update();
	TPaletteAxis *palette;
	palette = (TPaletteAxis*)histo->GetListOfFunctions()->FindObject("palette");
	// somehow, the palette->SetX2 is not working
	palette->SetX2NDC((412-gPad->GetX1())/(gPad->GetX2()-gPad->GetX1()));

	double ratio;	// must be 400/50 = 8
	ratio = fabs((double)(gPad->XtoPixel(400)-gPad->XtoPixel(0))/(gPad->YtoPixel(50)-gPad->YtoPixel(0)));

	c2->SetWindowSize(c2->GetWw()*8.0/ratio+offset_width, c2->GetWh()+offset_height); 
	c2->SetFixedAspectRatio();

	gPad->Modified();
	gPad->Update();
 	
 	
 	config.tfile->cd();
	char* fullName = new char[500];
	if (1/*organizeOutput*/) {
		config.tfile->mkdir(this->name);
		config.tfile->cd(this->name);
		sprintf(fullName, "%s", histoname);
	} else {
		sprintf(fullName, "%s-%s", this->name, histoname);
	}
	for (int ii = 0; ii < 500; ii++) { // Histo not callable from the cint cli if '-' in the name
		if (fullName[ii] == '-') {
			fullName[ii] = '_';
		}
		if (fullName[ii] == '\0') {
			break;
		}
	}
	c2->SetName(fullName);
	c2->Write(fullName);
	char *fileName = new char[500];
	sprintf(fileName, "%s%s-%s-%s%s", config.outPath, config.name, this->name,
			histoname, ".pdf");

	c2->SaveAs(fileName, "QQ");
	
	delete c2;
 	delete[] fileName;
	delete[] fullName;
	
 	gROOT->SetBatch(groot_batchmode);
 	gStyle->SetOptStat(gstyle_opt_stat);
	
	return;
}

/*
// geometrical cut for FEI3 VTT modules of TBNov2013
bool Botho::check_VTT_FEI3_TBNov2013(const Event &event) {
	// check if module is FEI3
	if(isFEI3) {
		// cut outermost columns (they are 600um instead of 400um)
		// -> pitchX/2 < event.trackX < (16+0.5)*pitchX
		// cut top 7 rows, because they are ganged with extra pixels outside
		// FE chip region
		// -> -pitchY/2 < event.trackY < (160-7)*pitchY
		if(event.trackX < 200 || event.trackX > 6600 ||
		   event.trackY < -25 || event.trackY > 7650) {
			return false;
		} else {
			return true;
	} else {	// if FEI4
		return check_fTrackCentralRegion(event);
	}
} */

// geometrical cut for FEI3 VTT modules of TBNov2013
event::eventflag_t Botho::check_central_region_VTT_FEI3_TBNov2013(const Event &event) {
	// cut outermost columns (they are 600um instead of 400um)
	// -> pitchX/2 < event.trackX < (16+0.5)*pitchX
	// cut top 7 rows, because they are ganged with extra pixels outside
	// FE chip region
	// -> -pitchY/2 < event.trackY < (160-7)*pitchY-pitchY/2
	if(event.trackX <= 200 || event.trackX >= 6600 ||
	   event.trackY <= -25 || event.trackY >= 7625) {
		return event::kBad;
	} else {
		return event::kGood;
	}
}

// second geometrical cut for FEI3 VTT modules of TBNov2013
// this cuts so that at least one regular row is outside of the cut area
// also in the current reconstruction the 0. column and 0. row are missing
// this means one row on the top and two on the bottom
// and two columns on every side
event::eventflag_t Botho::check_central_region_VTT_FEI3_TBNov2013_2(const Event &event) {
	// cut outermost two columns       outermost are irregular (they are 600um instead of 400um)
	// -> (2-0.5)*pitchX < event.trackX < (16-0.5)*pitchX
	// cut top 7 rows, because they are ganged with extra pixels outside
	// FE chip region
	// -> (2-0.5)*pitchY < event.trackY < (159-7-1-0.5)*pitchY
	if(event.trackX <= 600 || event.trackX >= 6200 ||
	   event.trackY <= 75 || event.trackY >= 7525) {
		return event::kBad;
	} else {
		return event::kGood;
	}
}

// geometrical and mask cut for FEI4 CIS reference module for GEOID100 of TBNov2013
event::eventflag_t Botho::check_fTrackRegion_geoid100_referenceDUT(const Event &event) {
	// cut all masked pixels and tracks in
	// columns and rows outside of beam area
	// triggered by the other DUT with a safety margin
	
/*	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
	                       // of matched cluster in vector of clusters of event
	if (matched == -1) {  // in case there is no matched cluster
		return event.fTrackRegion;	// leave event
	}
		
	// get reference to the matched cluster
	const vector<PllHit*>& matched_cluster = event.clusters.at(matched);
*/
	const double pitchX = 250;
	const double pitchY = 50;
	// good area col 15 to col 40 inclusively
	// and row 62 through 209 inclusively
	if(event.trackX < (pitchX*15)-pitchX/2 || event.trackX > (pitchX*40)+pitchX/2 ||
	   event.trackY < (pitchY*62)-pitchY/2 || event.trackY > (pitchY*209)+pitchY/2) {
		return event::kBad;
	} else {
		return event.fTrackRegion;
	}
	
//	return event::kGood;	
}

// geometrical and mask cut for VTT 200um FEI4 module (NP2-20-E4) for GEOIDs 198 to 204 of TBFeb2014
event::eventflag_t Botho::check_fTrackRegion_geoid198to204_DUT(const Event &event) {
	// cut all masked pixels and tracks in
	// columns and rows outside of beam area
	// triggered by the other DUT with a safety margin

	const double pitchX = 250;
	const double pitchY = 50;
	// good area col 21 to col 78 inclusively
	// and row 91 through 264 inclusively
	if(event.trackX < (pitchX*21)-pitchX/2 || event.trackX > (pitchX*78)+pitchX/2 ||
	   event.trackY < (pitchY*91)-pitchY/2 || event.trackY > (pitchY*264)+pitchY/2) {
		return event::kBad;
	} else {
		return event.fTrackRegion;
	}
	
//	return event::kGood;	
}