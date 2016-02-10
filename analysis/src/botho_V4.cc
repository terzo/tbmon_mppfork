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
	
	///////////////////////////
	// INITIALIZE HISTOGRAMS -------------------------------------------------
	/// define and initialize histograms
	
	/// cluster size histograms
	h_clusize_all = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
	h_clusize_matched = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
	h_clusize_unmatched = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
	h_clusize_matched_good_region = new TH1D("", ";cluster size;entries",
			n_rows, 0, n_rows);
	
	/// cluster multiplicity histograms
	h_cluster_multiplicity = new TH1D(
			"", ";number of clusters;events", 10, 0, 10);
	h_cluster_multiplicity_matched = new TH1D(
			"", ";number of clusters;events",10, 0, 10);
			
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
			
	/// cluster track maps
	for(int i=0; i<num_of_cs_track_histos; i++)
		h_track_cs[i] = new TH2D("", ";x [#mum];y [#mum]",
			(n_cols+4)*8, -2*pitchX, (n_cols+2)*pitchX,
			(n_rows+32) , -16*pitchY, (n_rows+16)*pitchY);

	/// cluster track maps projected into single pixel
	for(int i=0; i<num_of_cs_track_histos; i++)
		h_pixel_track_cs[i] = new TH2D("", ";x [#mum];y [#mum]",
			80, 0, pitchX,
			10, 0, pitchY);
	
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
		h_rawHits->Fill(col, row);
		int lv1 = (*it)->lv1;
		// fill h_lv1excluded if raw hit is not excluded by mask, but by lv1
		if((h_masks->GetCellContent(col, row) == 0)
				&& ((lv1 < lv1Min) || (lv1 > lv1Max))) {
			h_lv1excluded->Fill(col, row);
		}
	}
	
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
	
	/*// hit map
	for(vector<PllHit*>::const_iterator it = event.hits.begin();
			it != event.hits.end(); it++) {
		h_hitmap_raw_pos->Fill(event.trackX, event.trackY);
	} */
	

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
	
	// cluster multiplicity histogram
	h_cluster_multiplicity->Fill(event.clusters.size());	// increase entry of
                                        // corresponding number of clusters by 1
	
	// cluster size all histogram
	for(vector<event::cluster_t>::const_iterator it = event.clusters.begin(); 
	          it != event.clusters.end(); it++)	{ // loop through all clusters 
	                                               // in event
		h_clusize_all->Fill(it->size());
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
	   
	
	//////////////////////////////////
	/// CUT 2: GET MATCHED CLUSTER ///
	//////////////////////////////////
	
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);	// get number
	                       // of matched cluster in vector of clusters of event
	                       
	if (matched == -1) {  // in case there is no matched cluster
		// fill unmatched cluster histogram
		for(vector<event::cluster_t>::const_iterator
		     it = event.clusters.begin(); it != event.clusters.end(); it++) {
			h_clusize_unmatched->Fill(it->size());	// fill unmatched histogram
		}
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
	
	// matched cluster multiplicity histogram
	h_cluster_multiplicity_matched->Fill(event.clusters.size());
	
	// unmatched cluster size histogram
	for(vector<event::cluster_t>::const_iterator it = event.clusters.begin();
	     it != event.clusters.end(); it++) {	// loop over all clusters in event
		h_clusize_unmatched->Fill(it->size());	// unmatched histogram (matched
		                                   // ones will be removed after loop)
	}
	// fill histogram cluster size matched
	h_clusize_matched->Fill(matched_cluster.size());
	// decrease unmatched histogram by value of matched cluster
	h_clusize_unmatched->Fill(matched_cluster.size(),-1.0);
	
	// hit multiplicity histogram matched
	h_hit_multiplicity_matched->Fill(event.hits.size());
	
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
			
	
	// THIS PART WAS COPIED FROM RESIDUALS ANALYSIS BUT IS REMOVED HERE
/*	// Only use clusters of minimum total ToT min_sum_of_tot_per_cluster
	int min_sum_of_tot_per_cluster = 3;
	if (cluster::getSumTot(event.clusters.at(matched)) < min_sum_of_tot_per_cluster) {
		TBALOG(kDEBUG2) << "Sum of ToT entries in cluster too low: " << cluster::getSumTot(event.clusters.at(matched)) << " < " << min_sum_of_tot_per_cluster << "." << endl;
		return;
	}	*/
	
	
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
	
	
	// residuals
	double posX = cluster::getX(matched_cluster, event, cluster::kUnweighted);
	double posY = cluster::getY(matched_cluster, event, cluster::kUnweighted);
	  
	if (!(posX<0 || posY<0)){
		h_resX->Fill(event.trackX - posX);
		h_resY->Fill(event.trackY - posY);
	}
}


////////////////////////////////////////////////////////////////////////////////
/// FINALIZE ANALYSIS AFTER COMPLETE PROCESSING ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Botho::finalize(const TbConfig &config) {

	/////////////////////////////////////////
	/// set viewing ranges for histograms ///
	/////////////////////////////////////////
	h_clusize_all->GetXaxis()->SetRange(1, 6);
	h_clusize_matched->GetXaxis()->SetRange(1, 6);
	h_clusize_unmatched->GetXaxis()->SetRange(1, 6);
	h_clusize_matched_good_region->GetXaxis()->SetRange(1, 6);
	
	h_hit_multiplicity->GetXaxis()->SetRange(0, 10);
	h_hit_multiplicity_matched->GetXaxis()->SetRange(0, 10);


	/////////////////////////////////////////
	/// print and save histograms ///////////
	/////////////////////////////////////////
	
	// cluster histograms
	config.drawAndSave(this->name, "clusize_all", h_clusize_all);
	config.drawAndSave(this->name, "clusize_matched", h_clusize_matched);
	config.drawAndSave(this->name, "clusize_unmatched", h_clusize_unmatched);
	config.drawAndSave(this->name, "clusize_matched_good_region",
	          h_clusize_matched_good_region);
	
	// cluster multiplicities
	config.drawAndSave(this->name, "cluster_multiplicity",
	          h_cluster_multiplicity);
	config.drawAndSave(this->name, "cluster_multiplicity_matched",
	          h_cluster_multiplicity_matched);
	          
	// idens of hits (just for checking, whether all the hits in one analysis
	// actually have the same DUT iden (it is the DUT number specified in driver.cc))
	config.drawAndSave(this->name, "hit_iden", h_hit_iden);
	
	// hit multiplicities
	config.drawAndSave(this->name, "hit_multiplicity",
			h_hit_multiplicity);
	config.drawAndSave(this->name, "hit_multiplicity_matched",
			h_hit_multiplicity_matched);
	
	// hitmaps
	config.drawAndSave(this->name, "rawHits", h_rawHits);
	config.drawAndSave(this->name, "hitmap_all", h_hitmap_all);
	config.drawAndSave(this->name, "hitmap_all_check", h_hitmap_all_check);
	config.drawAndSave(this->name, "hitmap_raw_pos", h_hitmap_raw_pos);
	config.drawAndSave(this->name, "hitmap_all_pos", h_hitmap_all_pos);
	
	// outside lvl1cut
	config.drawAndSave(this->name, "lvl1_excluded", h_lv1excluded);
	
	// mask map
	config.drawAndSave(this->name, "masks", h_masks);
	
	// beam profile (different plots could be done with rebinning
	config.drawAndSave(this->name, "beamProfile10", h_beamProfile10);
	config.drawAndSave(this->name, "beamProfile5", h_beamProfile5);
	config.drawAndSave(this->name, "beamProfile2", h_beamProfile2);
	
	// angle profile
	h_angleProfile->Divide(h_angleProfile_numEvents);	// average angle profile
	config.drawAndSave(this->name, "angleProfile", h_angleProfile);
	
	// cluster track maps
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "track_cs", i+1);
		config.drawAndSave(this->name, histoname, h_track_cs[i]);
	}
	config.drawAndSave(this->name, "track_cs_overflow",
			h_track_cs[num_of_cs_track_histos-1]);
	
	// cluster track maps projected into single pixel
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d", "pixel_track_cs", i+1);
		config.drawAndSave(this->name, histoname, h_pixel_track_cs[i]);
	}
	config.drawAndSave(this->name, "pixel_track_cs_overflow",
			h_pixel_track_cs[num_of_cs_track_histos-1]);
			
	// rebin cluster track maps projected into single pixel
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		h_pixel_track_cs[i]->Rebin2D(2,2);
	}
	for(int i=0; i < num_of_cs_track_histos-1; i++) {
		char histoname[300];
		sprintf(histoname,"%s%d%s", "pixel_track_cs", i+1, "_rebin");
		config.drawAndSave(this->name, histoname, h_pixel_track_cs[i]);
	}
	config.drawAndSave(this->name, "pixel_track_cs_overflow_rebin",
			h_pixel_track_cs[num_of_cs_track_histos-1]);
 	
 	// residuals
 	saveToFile(config, this->name, "residuals", "resX", h_resX);
 	saveToFile(config, this->name, "residuals", "resY", h_resY);
 	
 	// finalize ToT analysis
 	finalize_tot_analysis(config);
	
	// finalize cuts analysis
	finalize_cut_analysis(config);
	
	// finalize efficiency analysis
	finalize_efficiency_analysis(config);

	
	////////////////////////////////
	/// free space of histograms ///
	////////////////////////////////
	
	// cluster histograms
	delete h_clusize_all;
	delete h_clusize_matched;
	delete h_clusize_unmatched;
	delete h_clusize_matched_good_region;
	
	delete h_cluster_multiplicity;
	delete h_cluster_multiplicity_matched;
	
	delete h_hit_iden;
	
	delete h_hit_multiplicity;
	delete h_hit_multiplicity_matched;

	delete h_rawHits;	
	delete h_hitmap_all;
	delete h_hitmap_all_check;
	delete h_hitmap_raw_pos;
	delete h_hitmap_all_pos;
	
	delete h_lv1excluded;
	delete h_masks;
	
	delete h_beamProfile10;
	delete h_beamProfile5;
	delete h_beamProfile2;
	
	delete h_angleProfile;
	delete h_angleProfile_numEvents;
	
	for(int i=0; i < num_of_cs_track_histos; i++)
		delete h_track_cs[i];
	for(int i=0; i < num_of_cs_track_histos; i++)
		delete h_pixel_track_cs[i];
	
	// residuals
	delete h_resX;
	delete h_resY;
}

////////////////////////////////////////////////////////////////////////////////
/// TOT - CHARGE ANALYSIS   ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// initialize tot analysis
void Botho::init_tot_analysis(const TbConfig &config) {
	// copy DUT information from the config object to members of this object
	const DUT* dut = config.getDut(this->iden);
	
	// get pitches
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();
	
	// get column and row numbers
	n_cols = dut->getNcols();
	n_rows = dut->getNrows();
	
	/// ToT module map
	h_module_tot_avg = new TH2D("", ";x;y;ToT", n_cols, 0, n_cols,
			n_rows, 0, n_rows);
	h_module_tot_hits = new TH2D("", ";x;y;hits", n_cols, 0, n_cols,
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
			
	/// residual histograms
	h_resX = new TH1D("", ";x [#mum]", 100, -pitchX * 2, pitchX * 2);
	h_resY = new TH1D("", ";y [#mum]", 100, -pitchY * 2, pitchY * 2);
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
	
	// ToT module histograms (unweighted algorithm for cluster position)
	h_module_tot_avg->Fill(
			cluster::getUnWeightedCol(matched_cluster),
			cluster::getUnWeightedRow(matched_cluster),
			cluster::getSumTot(matched_cluster));
	h_module_tot_hits->Fill(
			cluster::getUnWeightedCol(matched_cluster),
			cluster::getUnWeightedRow(matched_cluster));
			
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
	// ToT module maps
	saveToFile(config, this->name, "tot_analysis", "module_tot_sum",
			h_module_tot_avg);
	h_module_tot_avg->Divide(h_module_tot_hits);	// calculate average ToT
	saveToFile(config, this->name, "tot_analysis", "module_tot_hits",
			h_module_tot_hits);
	saveToFile(config, this->name, "tot_analysis", "module_tot_avg",
			h_module_tot_avg);
	
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
void Botho::init_cut_analysis(const TbConfig &config) {
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
void Botho::init_efficiency_analysis(const TbConfig &config) {
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
			
			
	/// efficiency.cc analysis rebuild ----------------------------------------
	effcc_ntracks = 0;
	effcc_nhits = 0;
	effcc_anyhit = 0;
	
	h_effcc_trackMap = new TH2D("", "Track Map;Column;Row", n_cols, 0, n_cols, n_rows, 0,
			n_rows);
	h_effcc_hitMap = new TH2D("", "Hit Map;Column;Row", n_cols, 0, n_cols, n_rows, 0,
			n_rows);
	h_effcc_anyHitMap = new TH2D("", "Any Hit Map;Column;Row", n_cols, 0, n_cols, n_rows,
			0, n_rows);

}			

// analyze
void Botho::efficiency_analysis(const TbConfig &config, const Event &event) {
	
	////////////////////////////////////////////
	/// RAW DATA ANALYSIS                    ///
	////////////////////////////////////////////
	
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
	
	/// TODO
	/// redo analysis as in efficiency.cc
	
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
 	double rebinX = 10;
 	double rebinY = 5;
 	
 	h_trackPixelMap->Rebin2D(rebinX,rebinY);
 	h_hitPixelMap->Rebin2D(rebinX,rebinY);
 	
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
}	

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

