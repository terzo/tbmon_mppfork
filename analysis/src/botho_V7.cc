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
	h_profileX = new TH1D("", ";x [#mum]; hits",
			n_cols*pitchX/10, 0, n_cols*pitchX);
	h_profileY = new TH1D("", ";y [#mum]; hits",
			n_rows*pitchY/10, 0, n_rows*pitchY);
	h_profileX_pixel = new TH1D("", ";x [#mum]; hits",
			pitchX, 0, pitchX);
	h_profileY_pixel = new TH1D("", ";y [#mum]; hits",
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
			
	/// beam profiles	
	res1 = 10;
	h_beamProfile10 = new TH2D("", ";x [um];y [um]",
			(n_cols*3)*res1, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1, -1*n_rows*pitchY, (n_rows*2)*pitchY);
	res1 = 5;
	h_beamProfile5 = new TH2D("", ";x [um];y [um]",
			(n_cols*3)*res1, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1, -1*n_rows*pitchY, (n_rows*2)*pitchY);
	res1 = 2;
	h_beamProfile2 = new TH2D("", ";x [um];y [um]",
			(n_cols*3)*res1, -1*n_cols*pitchX, (n_cols*2)*pitchX,
			(n_rows*3)*res1, -1*n_rows*pitchY, (n_rows*2)*pitchY);
			
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
	h_profileX_pixel->Fill(tbutils::getFoldedX(event));
	h_profileY_pixel->Fill(tbutils::getPixelY(event));
	
	// hit map
	h_hitmap_raw_pos->Fill(event.trackX, event.trackY);
	
	// beam profiles
	h_beamProfile10->Fill(event.trackX, event.trackY);
	h_beamProfile5->Fill(event.trackX, event.trackY);
	h_beamProfile2->Fill(event.trackX, event.trackY);
	
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
	
	
	/// CUT 3:  EVENT.FTRACKREGION == EVENT::KGOOD ///
	
	// Track in good region? (central and unmasked) (see event.h)
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	
	// residuals
	double posX = cluster::getX(matched_cluster, event, cluster::kUnweighted);
	double posY = cluster::getY(matched_cluster, event, cluster::kUnweighted);
	  
	if (!(posX<0 || posY<0)){
		h_resX->Fill(event.trackX - posX);
		h_resY->Fill(event.trackY - posY);
	}
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
	
	// hit profiles
	saveToFile(config, this->name, "hit_analysis", "profileX", h_profileX);
	saveToFile(config, this->name, "hit_analysis", "profileY", h_profileY);
	saveToFile(config, this->name, "hit_analysis", "profileX_pixel",
			h_profileX_pixel);
	saveToFile(config, this->name, "hit_analysis", "profileY_pixel",
			h_profileY_pixel);
	
	// beam profile (different plots could be done with rebinning
	saveToFile(config, this->name, "hit_analysis", "beamProfile10", h_beamProfile10);
	saveToFile(config, this->name, "hit_analysis", "beamProfile5", h_beamProfile5);
	saveToFile(config, this->name, "hit_analysis", "beamProfile2", h_beamProfile2);
	
	// angle profile
	h_angleProfile->Divide(h_angleProfile_numEvents);	// average angle profile
	saveToFile(config, this->name, "hit_analysis", "angleProfile", h_angleProfile);
 	
 	// residuals
 	saveToFile(config, this->name, "hit_analysis", "resX", h_resX);
 	saveToFile(config, this->name, "hit_analysis", "resY", h_resY);
 	
 	
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
	
	delete h_beamProfile10;
	delete h_beamProfile5;
	delete h_beamProfile2;
	
	delete h_angleProfile;
	delete h_angleProfile_numEvents;
	
	// residuals
	delete h_resX;
	delete h_resY;
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
			n_rows, 0, n_rows);
	h_clusize_matched = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
	h_clusize_unmatched = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
	h_clusize_matched_good_region = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
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
	
	// cluster track maps
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
	
	/// CUT 3:  EVENT.FTRACKREGION == EVENT::KGOOD ///
	// Track in good region? (central and unmasked) (see event.h)
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	
	// fill diagram cluster size matched and good region
	h_clusize_matched_good_region->Fill(event.clusters.at(matched).size());
	
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
	if(matched_cluster.size() < num_of_cs_track_histos) {	// check if cluster size smaller than constant
		if(matched_cluster.size() > 0) {	// check whether cluster size positive
			// fill cluster size histograms
			h_pixel_track_cs[matched_cluster.size()-1]->Fill(trackX, trackY);
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
	}
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
	
	// cluster track maps projected into single pixel
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "pixel_track_cs", i+1);
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs[i]);
	}
	saveToFile(config, this->name, "cluster_analysis", "pixel_track_cs_overflow",
			h_pixel_track_cs[num_of_cs_track_histos-1]);
			
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		// rebin cluster track maps projected into single pixel
		h_pixel_track_cs[i]->Rebin2D(2,2);
		char histoname[300];
		sprintf(histoname,"%s%d%s", "pixel_track_cs", i+1, "_rebin");
		saveToFile(config, this->name, "cluster_analysis", histoname,
				h_pixel_track_cs[i]);
	}
	h_pixel_track_cs[num_of_cs_track_histos-1]->Rebin2D(2,2);
	saveToFile(config, this->name, "cluster_analysis",
			"pixel_track_cs_overflow_rebin",
			h_pixel_track_cs[num_of_cs_track_histos-1]);
			
	// delete cluster histograms
	delete h_clusize_all;
	delete h_clusize_matched;
	delete h_clusize_unmatched;
	delete h_clusize_matched_good_region;
	
	delete h_cluster_multiplicity;
	delete h_cluster_multiplicity_matched;
	
	for(int i=0; i < num_of_cs_track_histos; i++)
		delete h_track_cs[i];
	for(int i=0; i < num_of_cs_track_histos; i++)
		delete h_pixel_track_cs[i];
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
	} else if(n_cols==18 && n_rows==160) {	// case of FE-I3
		n_tot_bins = 256;
	} else {
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
	  (event.fTrackRegion == event::kGood) &&
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

