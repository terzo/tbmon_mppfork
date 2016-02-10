#include "residuals.h"

void Residuals::init(TbConfig &config) {

	// copy information from the config object to member of this object
	const DUT* dut = config.getDut(this->iden);
	pitchX = dut->getPitchX();
	pitchY = dut->getPitchY();

	// get options or set standard values
	use_minimum_sum_of_tot_per_cluster =
			(bool) config.cmdLineExtras_argGetter(
					(string) "A_residuals_UseMinSumTotClus", (bool) true,
					(string) "Use cut: Minimum sum of a cluster to be taken into account for residual plots.");
	minimum_sum_of_tot_per_cluster =
			(int) config.cmdLineExtras_argGetter(
					(string) "A_residuals_MinSumTotClus", (double) 3,
					(string) "Minimum sum of a cluster to be taken into account for residual plots.");
	use_minimum_entries_to_plot_histo =
			(bool) config.cmdLineExtras_argGetter(
					(string) "A_residuals_UseMinEntrToPlotHisto", (bool) false,
					(string) "Use cut: Minimum number of entries in a residual histogram, so that it is actually plotted.");
	minimum_entries_to_plot_histo =
			(int) config.cmdLineExtras_argGetter(
					(string) "A_residuals_MinEntrToPlotHisto", (double) 20,
					(string) "Minimum number of entries in a residual histogram, so that it is actually plotted.");

	// initialize standard histograms (all other histograms are created dynamically as need arises)
	h_angleX = new TH1D("", ";Beam dx/dz;Entries", 1000, -0.001, 0.001);
	h_angleY = new TH1D("", ";Beam dy/dz;Entries", 1000, -0.001, 0.001);
	h_clusizeX = new TH1D("", ";Cluster span [col];Entries", dut->getNcols(), 0,
			dut->getNcols());
	h_clusizeY = new TH1D("", ";Cluster span [row];Entries", dut->getNrows(), 0,
			dut->getNrows());
	h_clusize = new TH1D("", ";Cluster size;Entries", dut->getNrows(), 0,
			dut->getNrows());
	h_clusizeXY = new TH2D("", ";Cluster span [col];Cluster span [row]",
			dut->getNcols(), 0, dut->getNcols(), dut->getNrows(), 0,
			dut->getNrows());
	h_ecorrsX = new TH1D("", ";Eta Correction X [#mum];Entries", 100, 0,
			pitchX);
	h_ecorrsY = new TH1D("", ";Eta Correction Y [#mum];Entries", 100, 0,
			pitchY);

}

void Residuals::getResiduals(const vector<PllHit*> &cluster, const Event &event,
		vector<TH1D*> &histosX, vector<TH1D*> &histosY,
		vector<TH1D*> &histosAll) {
      
	for (int alg=0; alg<cluster::kUnknown; alg++){
	  double posX = cluster::getX(cluster, event, alg);
	  double posY = cluster::getY(cluster, event, alg);
	  
	  if (posX<0 || posY<0){
	    continue;
	  }
	  histosX[alg]->Fill(event.trackX - posX);
	  histosY[alg]->Fill(event.trackY - posY);
	  
	  double distance=sqrt(pow(event.trackX - posX, 2) + pow(event.trackY - posY, 2));
	  histosAll[alg]->Fill(distance);
	} 
}

void Residuals::event(const TbConfig &config, const Event &event) {
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		return;
	}
	TBALOG(kDEBUG2) <<
			"Event passed first cut: \"event.fTrack == event::kGood\"." << endl;

	h_angleX->Fill(event.dxdz);
	TBALOG(kDEBUG2) << "Filled dx/dz= " << event.dxdz << " into histogram."
			<< endl;

	h_angleY->Fill(event.dydz);
	TBALOG(kDEBUG2) << "Filled dy/dz= " << event.dydz << " into histogram."
			<< endl;

	// Are we in a good region of the chip?
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	TBALOG(kDEBUG2) <<
			"Event passed second cut: \"event.fTrackRegion == event::kGood\"."
			<< endl;

	//Check that we have one and only one cluster
	//Why should this be necessary? Just take the cluster which works best with the matching criterium and be done.
	//if( event.clusters.size() != 1) { return; }
	//TBALOG(kDEBUG2) << "Event passed third cut: \"event.clusters.size() == 1\"." <<  endl;
	TBALOG(kDEBUG2) <<
			"Cut: \"event.clusters.size() == 1\" is no longer applied." << endl;

	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);
	TBALOG(kDEBUG2) << "Found " << matched << " matched clusters." << endl;
	if (matched == -1) {
		TBALOG(kDEBUG2) << "No matched clusters --> leaving event." << endl;
		return;
	}

	const vector<PllHit*>& cluster = event.clusters.at(matched);
	// Check that cluster has charge
	if (use_minimum_sum_of_tot_per_cluster == true) {
		if (cluster::getSumTot(cluster) < minimum_sum_of_tot_per_cluster) {
			TBALOG(kDEBUG2) << "Sum of ToT entries in cluster too low: "
					<< cluster::getSumTot(cluster) << " < "
					<< minimum_sum_of_tot_per_cluster << "." << endl;
			return;
		}; // make FE-dependent ToT cut?
	}

	// Plot residuals based on cluster size
	int nCols = cluster::spannedCols(cluster);
	int nRows = cluster::spannedRows(cluster);
	h_clusizeX->Fill(nCols);
	h_clusizeY->Fill(nRows);
	h_clusize->Fill(cluster.size());
	h_clusizeXY->Fill(nCols, nRows);

	// check whether there was already created a histogram for the cluster size nCols in X --> this cluster size already occurred and we can just fill the histogram
	if (hhh_resX[nCols].empty()) {
		for (int method = 0; method < cluster::kUnknown; method++) {
			char title[400];
			sprintf(title, "%i-cluSize%i-resX;X Residual [#mum];Entries", iden,
					nCols);
			hhh_resX[nCols].push_back(
					new TH1D("", title, 100, -pitchX * nCols, pitchX * nCols));
			TBALOG(kDEBUG) << "initiated " << title << "," << method << endl;
		}
	}
	// check whether there was already created a histogram for the cluster size nRows in Y --> this cluster size already occurred and we can just fill the histogram
	if (hhh_resY[nRows].empty()) {
		for (int method = 0; method < cluster::kUnknown; method++) {
			char title[400];
			sprintf(title, "%i-cluSize%i-resY;Y Residual [#mum];Entries", iden,
					nRows);
			hhh_resY[nRows].push_back(
					new TH1D("", title, 100, -(pitchY * 2 * nRows),
							pitchY * 2 * nRows));
			TBALOG(kDEBUG) << "initiated " << title << "," << method << endl;
		}
	}
	// check whether there was already created a histogram for the total cluster size  cluster.size() --> this cluster size already occurred and we can just fill the histogram
	if (hhh_res[cluster.size()].empty()) {
		for (int method = 0; method < cluster::kUnknown; method++) {
			char title[400];
			sprintf(title, "%i-cluSize%i-res;Residual [#mum];Entries", iden,
					(int) cluster.size());
			hhh_res[cluster.size()].push_back(
					new TH1D("", title, 100,
							-(max(2 * pitchY * nRows, pitchX * nCols)),
							max(2 * pitchY * nRows, pitchX * nCols)));
			TBALOG(kDEBUG) << "initiated " << title << "," << method << endl;
		}
	}

	getResiduals(cluster, event, hhh_resX[nCols], hhh_resY[nRows],
			hhh_res[(int) cluster.size()]);

	TBALOG(kDEBUG2) << "Filled general histograms." << endl;

	//WARNING: implemented eta correction can only handle 1 and 2 hit clusters --> only handle appropriate cluster sizes
	if (cluster.size() <= 2) {
		//Tight matching for numbers to quote
		if (fabs(event.trackY 
			  - cluster::getEtaCorrectedY(cluster, event, event.dut))
		      < pitchY) {
			if (hhh_resX[0].empty()) {
				for (int method = 0; method < cluster::kUnknown; method++) {
					char title[400];
					sprintf(title, "%i-cluSize0-resX;X Residual [#mum];Entries",
							iden);
					hhh_resX[0].push_back(
							new TH1D("", title, 100, -pitchX * nCols,
									pitchX * nCols));
					TBALOG(kDEBUG) << "initiated " << title << "," << method
							<< endl;
				};
			}
			if (hhh_resY[0].empty()) {
				for (int method = 0; method < cluster::kUnknown; method++) {
					char title[400];
					sprintf(title, "%i-cluSize0-resY;Y Residual [#mum];Entries",
							iden);
					hhh_resY[0].push_back(
							new TH1D("", title, 100, -(2 * pitchY * nRows),
									2 * pitchY * nRows));
					TBALOG(kDEBUG) << "initiated " << title << "," << method
							<< endl;
				};
			}
			if (hhh_res[0].empty()) {
				for (int method = 0; method < cluster::kUnknown; method++) {
					char title[400];
					sprintf(title, "%i-cluSize0-res;Residual [#mum];Entries",
							iden);
					hhh_res[0].push_back(
							new TH1D("", title, 100,
									-(max(2 * pitchY * nRows, pitchX * nCols)),
									max(2 * pitchY * nRows, pitchX * nCols)));
					TBALOG(kDEBUG) << "initiated " << title << "," << method
							<< endl;
				};
			}
			getResiduals(cluster, event, hhh_resX[0], hhh_resY[0], hhh_res[0]);
		}

		if (event.fEtaCorrections == event::kGood) {
			if (fabs(
					event.trackX
							- cluster::getEtaCorrectedX(cluster, event,
									event.dut)) > pitchX) {
				return;
			}
			if (fabs(
					event.trackY
							- cluster::getEtaCorrectedY(cluster, event,
									event.dut)) > pitchY) {
				return;
			}
			// All cluster sizes
			checkEtaCorrections(cluster, event);
		}
	}

}

void Residuals::checkEtaCorrections(const vector<PllHit*> &cluster,
		const Event &event) {
	double qmeanX = cluster::getChargeWeightedX(cluster, event);
	double qmeanY = cluster::getChargeWeightedY(cluster, event);
	double ecorrX = cluster::getEtaCorrectedX(cluster, event, event.dut);
	double ecorrY = cluster::getEtaCorrectedY(cluster, event, event.dut);
	if (ecorrX != qmeanX) {
		while (ecorrX > event.dut->getPitchX()) {
			ecorrX -= event.dut->getPitchX();
		}
		h_ecorrsX->Fill(ecorrX);
	}
	if (ecorrY != qmeanY) {
		while (ecorrY > event.dut->getPitchY()) {
			ecorrY -= event.dut->getPitchY();
		}
		h_ecorrsY->Fill(ecorrY);
	}

}

void Residuals::finalize(const TbConfig &config) {

	TBALOG(kINFO) << config.name << ":" << endl;
	config.drawAndSave(this->name, "angleX", h_angleX);
	config.drawAndSave(this->name, "angleY", h_angleY);
	config.drawAndSave(this->name, "clusizeX", h_clusizeX);
	config.drawAndSave(this->name, "clusizeY", h_clusizeY);
	config.drawAndSave(this->name, "clusize", h_clusize);
	config.drawAndSave(this->name, "clusizeXY", h_clusizeXY);
	config.drawAndSave(this->name, "etacorrsX", h_ecorrsX);
	config.drawAndSave(this->name, "etacorrsY", h_ecorrsY);

	map<int, vector<TH1D*> >::iterator it;
	int ii;

	for (it = hhh_resX.begin(); it != hhh_resX.end(); it++) {
		char* name = new char[500];
		ii = (*it).first;

		//write basic residual histogram in root file and separate files if desired
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			if ((use_minimum_entries_to_plot_histo == false)
					|| (hhh_resX[ii][jj]->GetEntries()
							>= minimum_entries_to_plot_histo)) {
				config.saveToFile(this->name,
						TString::Format(
								"residual-in-X-clusize-%i-with-%s-algorithm",
								ii, cluster::getAlgName(jj).c_str()).Data(),
						hhh_resX[ii][jj]);
				config.drawToFile(this->name,
						TString::Format(
								"residual-in-X-clusize-%i-with-%s-algorithm",
								ii, cluster::getAlgName(jj).c_str()).Data(),
						hhh_resX[ii][jj]);
			} else {
				TBALOG(kDEBUG) << "Residual plot for cluster size " << ii
						<< " and algorithm " << cluster::getAlgTitle(jj)
						<< " not drawn as there are too few entries ("
						<< hhh_resX[ii][jj]->GetEntries()
						<< ") and A_residuals_UseMinEntrToPlotHisto is true."
						<< endl;
			}
		}

		//dump basic information on residuals of that cluster size
		TBALOG(kINFO) << "X(cluSize " << ii << "):" << endl;
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			TBALOG(kINFO) << " " << cluster::getAlgTitle(jj) << ": "
					<< hhh_resX[ii][jj]->GetRMS() << " ("
					<< hhh_resX[ii][jj]->GetEntries() << " Entries)" << endl;
		}
		TBALOG(kINFO) << endl;

		//write comparison plots of residuals with different cluster center finding algorthm in root file and separate files if desired
		compareResidualsForOneClustersize(config, this->name,
				TString::Format(
						"residual-in-X-clusize-%i-comparison-between-algorithm",
						ii), hhh_resX[ii]);
		TBALOG(kDEBUG) << "Compare Residuals for cluster size in x = " << ii
				<< " done." << endl;

		if (hhh_resX.find(ii) != hhh_resX.end()) {
			TBALOG(kDEBUG) <<
					"Entered fitting subroutine first check - for cluster size "
					<< ii << " in x an entry exists." << endl;
			for (int jj = 0; jj < cluster::kUnknown; ++jj) {
				if (hhh_resX[ii][jj] != 0) {
					TBALOG(kDEBUG) <<
							"Entered fitting subroutine second check - for cluster size "
							<< ii << " in x and " << cluster::getAlgTitle(jj)
							<< " algorithm a histogram exists. Handing over to fitting methods."
							<< endl;
					if ((use_minimum_entries_to_plot_histo == false)
							|| (hhh_resX[ii][jj]->GetEntries()
									>= minimum_entries_to_plot_histo)) {
						fitResiduals(config, this->name,
								TString::Format(
										"fitResidual-in-X-clusize-%i-algorithm-%s",
										ii, cluster::getAlgName(jj).c_str()),
								hhh_resX[ii][jj], pitchX / 2.0, -pitchX * ii,
								pitchX * ii);
						fitAllClusResiduals(config, this->name,
								TString::Format(
										"fitResidual-multiG-in-X-clusize-%i-algorithm-%s",
										ii, cluster::getAlgName(jj).c_str()),
								hhh_resX[ii][jj], pitchX / 2.0, -pitchX * ii,
								pitchX * ii);
						TBALOG(kDEBUG) << "Fitting for cluster size " << ii
								<< " in x is done." << endl;
					} else {
						TBALOG(kDEBUG) << "Residual fit for cluster size "
								<< ii << " and algorithm " << cluster::getAlgTitle(jj)
								<< " not drawn as there are too few entries ("
								<< hhh_resX[ii][jj]->GetEntries()
								<< ") and A_residuals_UseMinEntrToPlotHisto is true."
								<< endl;
					}
				} else {
					TBALOG(kDEBUG) <<
							"Entered fitting subroutine second check - for cluster size "
							<< ii << " in x and " << cluster::getAlgTitle(jj)
							<< " algorithm no histogram exists: skip." << endl;
				}
			}
		} else {
			TBALOG(kDEBUG) <<
					"Entered fitting subroutine first check - for cluster size "
					<< ii << " in x no entry exists: skip." << endl;
		}

		delete[] name;
	}

	for (it = hhh_resY.begin(); it != hhh_resY.end(); it++) {
		char* name = new char[500];
		ii = (*it).first;

		//write basic residual histogram in root file and separate files if desired
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			if ((use_minimum_entries_to_plot_histo == false)
					|| (hhh_resY[ii][jj]->GetEntries()
							>= minimum_entries_to_plot_histo)) {
				config.saveToFile(this->name,
						TString::Format(
								"residual-in-Y-clusize-%i-with-%s-algorithm",
								ii, cluster::getAlgName(jj).c_str()).Data(),
						hhh_resY[ii][jj]);
			} else {
				TBALOG(kDEBUG) << "Residual plot for cluster size " << ii
						<< " and algorithm " << cluster::getAlgTitle(jj)
						<< " not drawn as there are too few entries ("
						<< hhh_resY[ii][jj]->GetEntries()
						<< ") and A_residuals_UseMinEntrToPlotHisto is true."
						<< endl;
			}
		}

		//dump basic information on residuals of that cluster size
		TBALOG(kINFO) << "Y(cluSize " << ii << "):" << endl;
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			TBALOG(kINFO) << " " << cluster::getAlgTitle(jj) << ": "
					<< hhh_resY[ii][jj]->GetRMS() << " ("
					<< hhh_resY[ii][jj]->GetEntries() << " Entries)" << endl;
		}
		TBALOG(kINFO) << endl;

		//write comparison plots of residuals with different cluster center finding algorthm in root file and separate files if desired
		compareResidualsForOneClustersize(config, this->name,
				TString::Format(
						"residual-in-Y-clusize-%i-comparison-between-algorithm",
						ii), hhh_resY[ii]);
		TBALOG(kDEBUG) << "Compare Residuals for cluster size in y = " << ii
				<< " done." << endl;

		if (hhh_resY.find(ii) != hhh_resY.end()) {
			TBALOG(kDEBUG) <<
					"Entered fitting subroutine first check - for cluster size "
					<< ii << " in y an entry exists." << endl;
			for (int jj = 0; jj < cluster::kUnknown; ++jj) {
				if (hhh_resY[ii][jj] != 0) {
					TBALOG(kDEBUG) <<
							"Entered fitting subroutine second check - for cluster size "
							<< ii << " in y and " << cluster::getAlgTitle(jj)
							<< " algorithm a histogram exists. Handing over to fitting methods."
							<< endl;
					if ((use_minimum_entries_to_plot_histo == false)
							|| (hhh_resY[ii][jj]->GetEntries()
									>= minimum_entries_to_plot_histo)) {
						fitResiduals(config, this->name,
								TString::Format(
										"fitResidual-in-Y-clusize-%i-algorithm-%s",
										ii, cluster::getAlgName(jj).c_str()),
								hhh_resY[ii][jj], pitchY / 2.0, -pitchY * ii,
								pitchY * ii);
						fitAllClusResiduals(config, this->name,
								TString::Format(
										"fitResidual-multiG-in-Y-clusize-%i-algorithm-%s",
										ii, cluster::getAlgName(jj).c_str()),
								hhh_resY[ii][jj], pitchY / 2.0, -pitchY * ii,
								pitchY * ii);
						TBALOG(kDEBUG) << "Fitting for cluster size " << ii
								<< " in y is done." << endl;
					} else {
						TBALOG(kDEBUG) << "Residual fit for cluster size "
								<< ii << " and algorithm " << cluster::getAlgTitle(jj)
								<< " not drawn as there are too few entries ("
								<< hhh_resY[ii][jj]->GetEntries()
								<< ") and A_residuals_UseMinEntrToPlotHisto is true."
								<< endl;
					}
				} else {
					TBALOG(kDEBUG) <<
							"Entered fitting subroutine second check - for cluster size "
							<< ii << " in y and " << cluster::getAlgTitle(jj)
							<< " algorithm no histogram exists: skip." << endl;
				}
			}
		} else {
			TBALOG(kDEBUG) <<
					"Entered fitting subroutine first check - for cluster size "
					<< ii << " in y no entry exists: skip." << endl;
		}
		delete[] name;
	}

	for (it = hhh_res.begin(); it != hhh_res.end(); it++) {
		char* name = new char[500];
		ii = (*it).first;

		//write basic residual histogram in root file and separate files if desired
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			if ((use_minimum_entries_to_plot_histo == false)
					|| (hhh_res[ii][jj]->GetEntries()
							>= minimum_entries_to_plot_histo)) {
				config.saveToFile(this->name,
						TString::Format(
								"residual-in-total-clusize-%i-with-%s-algorithm",
								ii, cluster::getAlgName(jj).c_str()).Data(),
						hhh_res[ii][jj]);
			} else {
				TBALOG(kDEBUG) << "Residual plot for cluster size " << ii
						<< " and algorithm " << cluster::getAlgTitle(jj)
						<< " not drawn as there are too few entries ("
						<< hhh_res[ii][jj]->GetEntries()
						<< ") and A_residuals_UseMinEntrToPlotHisto is true."
						<< endl;
			}
		}

		//dump basic information on residuals of that cluster size
		TBALOG(kINFO) << "Total(cluSize " << ii << "):";
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			TBALOG(kINFO) << " " << cluster::getAlgTitle(jj) << ": "
					<< hhh_res[ii][jj]->GetRMS() << " ("
					<< hhh_res[ii][jj]->GetEntries() << " Entries)" << endl;
		}
		TBALOG(kINFO) << endl;

		//write comparison plots of residuals with different cluster center finding algorthm in root file and separate files if desired
		compareResidualsForOneClustersize(config, this->name,
				TString::Format(
						"residual-in-total-clusize-%i-comparison-between-algorithm",
						ii), hhh_res[ii]);
		TBALOG(kDEBUG) << "Compare Residuals for total cluster size = " << ii
				<< " done." << endl;

		if (hhh_res.find(ii) != hhh_res.end()) {
			TBALOG(kDEBUG) <<
					"Entered fitting subroutine first check - for total cluster size "
					<< ii << " an entry exists." << endl;
			for (int jj = 0; jj < cluster::kUnknown; ++jj) {
				if (hhh_res[ii][jj] != 0) {
					TBALOG(kDEBUG) <<
							"Entered fitting subroutine second check - for total cluster size "
							<< ii << " and " << cluster::getAlgTitle(jj)
							<< " algorithm a histogram exists. Handing over to fitting methods."
							<< endl;
					if ((use_minimum_entries_to_plot_histo == false)
							|| (hhh_res[ii][jj]->GetEntries()
									>= minimum_entries_to_plot_histo)) {
						fitResiduals(config, this->name,
								TString::Format(
										"fitResidual-in-total-clusize-%i-algorithm-%s",
										ii, cluster::getAlgName(jj).c_str()),
								hhh_res[ii][jj],
								max(pitchX / 2.0, pitchY / 2.0),
								-(max(pitchX * ii, 2 * pitchY * ii)),
								max(pitchX * ii, 2 * pitchY * ii));
						fitAllClusResiduals(config, this->name,
								TString::Format(
										"fitResidual-multiG-in-total-clusize-%i-algorithm-%s",
										ii, cluster::getAlgName(jj).c_str()),
								hhh_res[ii][jj],
								max(pitchX / 2.0, pitchY / 2.0),
								-(max(pitchX * ii, 2 * pitchY * ii)),
								max(pitchX * ii, 2 * pitchY * ii));
						TBALOG(kDEBUG) << "Fitting for total cluster size "
								<< ii << " is done." << endl;
					} else {
						TBALOG(kDEBUG) << "Residual fit for cluster size "
								<< ii << " and algorithm " << cluster::getAlgTitle(jj)
								<< " not drawn as there are too few entries ("
								<< hhh_res[ii][jj]->GetEntries()
								<< ") and A_residuals_UseMinEntrToPlotHisto is true."
								<< endl;
					}
				} else {
					TBALOG(kDEBUG) <<
							"Entered fitting subroutine second check - for total cluster size "
							<< ii << " and " << cluster::getAlgTitle(jj)
							<< " algorithm no histogram exists: skip." << endl;
				}
			}
		} else {
			TBALOG(kDEBUG) <<
					"Entered fitting subroutine first check - for total cluster size "
					<< ii << " no entry exists: skip." << endl;
		}
		delete[] name;
	}

	//clean up heap
	delete h_angleX;
	delete h_angleY;
	delete h_clusizeX;
	delete h_clusizeY;
	delete h_clusize;
	delete h_clusizeXY;
	delete h_ecorrsX;
	delete h_ecorrsY;

	for (it = hhh_resX.begin(); it != hhh_resX.end(); it++) {
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			delete hhh_resX[(*it).first][jj];
		}
	}

	for (it = hhh_resY.begin(); it != hhh_resY.end(); it++) {
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			delete hhh_resY[(*it).first][jj];
		}
	}

	for (it = hhh_res.begin(); it != hhh_res.end(); it++) {
		for (int jj = 0; jj < cluster::kUnknown; ++jj) {
			delete hhh_res[(*it).first][jj];
		}
	}
}

void Residuals::compareResidualsForOneClustersize(const TbConfig& config,
		const char* analysisName, const char* histoName, vector<TH1D*>& hh) {
	//Compare the different algorithms
	char* fileName = config.buildHistName(analysisName, histoName);
	char newHistoName[150];
	sprintf(newHistoName, "%s-%s", analysisName, histoName);

	
	double maxval = 0;
	for (int jj=0; jj<cluster::kUnknown; jj++){
	  if (hh[jj]->GetMaximum() > maxval){
	    maxval = hh[jj]->GetMaximum();
	  }
	}
	  
	double lw = 2.0;
	TCanvas* canvas = new TCanvas(TString::Format("can_%s", fileName),
			TString::Format("can_%s", fileName), 10, 20, 700, 500);
	
	bool firstPlot = true;
	for (int jj=0; jj<cluster::kUnknown; jj++){
	  hh[jj]->SetMaximum(maxval * 1.1);
	  hh[jj]->SetStats(false);
	  hh[jj]->SetLineWidth(lw);
	  hh[jj]->SetLineColor(cluster::lineColor[jj]);
	  hh[jj]->SetFillColor(cluster::areaColor[jj]);
	  if (cluster::areaColor[jj] != kWhite){
	    hh[jj]->SetFillStyle(3002);
	  }else{
	    hh[jj]->SetFillStyle(0);
	  }
	  if ((use_minimum_entries_to_plot_histo == false)
			|| (hh[jj]->GetEntries() >= minimum_entries_to_plot_histo)) {
	    if (firstPlot){
	      hh[jj]->Draw();
	      firstPlot = false;
	    }else{
	      hh[jj]->Draw("same");
	    }
	} else {
		TBALOG(kDEBUG) << "Residual plot for algorithm " << cluster::getAlgTitle(jj)
				<< " not drawn in comparison plot " << histoName
				<< " for this cluster size as there are too few entries and A_residuals_UseMinEntrToPlotHisto is true."
				<< endl;
	}
      }
	  
	hh[cluster::kUnweighted]->SetStats(true);


	TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89);
	leg->SetFillColor(kWhite);
	leg->SetBorderSize(1);
	
	for(int jj=0; jj<cluster::kUnknown; jj++){
	  if ((use_minimum_entries_to_plot_histo == false)
		|| (hh[jj]->GetEntries() >= minimum_entries_to_plot_histo)) {
	    if (cluster::areaColor[jj] != kWhite){
	      leg->AddEntry(hh[jj], cluster::getAlgTitle(jj).c_str(), "F");
	    }else{
	      leg->AddEntry(hh[jj], cluster::getAlgTitle(jj).c_str(), "L");
	    }  
	  }
	}
	leg->Draw();
	
	bool savePlot = (!use_minimum_entries_to_plot_histo);
	for (int jj=0; jj<cluster::kUnknown; jj++){
	  savePlot |= (hh[jj]->GetEntries() >= minimum_entries_to_plot_histo);
	}
	if (savePlot) {
		if (config.plotExtension != "none") {
			cout << "test" << endl;
			canvas->SaveAs(Form("%s/%s", config.outPath, newHistoName));
		}
		canvas->Write(newHistoName);
	} else {
		TBALOG(kDEBUG) << "Residual comparison plot " << histoName
				<< "  not drawn for this cluster size as there are too few entries in all plots and A_residuals_UseMinEntrToPlotHisto is true."
				<< endl;
	}

	delete[] fileName;
	delete canvas;
	delete leg;

	return;
}

void Residuals::fitResiduals(const TbConfig& config, const char* analysisName,
		const char* histoName, TH1 *h_residuals, double halfPitch, double xmin,
		double xmax) {
	char* fileName = config.buildHistName(analysisName, histoName);
	char newHistoName[150];
	sprintf(newHistoName, "%s-%s", analysisName, histoName);

	TF1 *func = new TF1("func", tbutils::fitFunc, xmin, xmax, 5);
	func->SetParameters(8.0, h_residuals->Integral("width"), 0.0, halfPitch,
			0.0);
	func->FixParameter(3, halfPitch);
	func->SetLineColor(2);
	h_residuals->Fit(func, "+ Q");

	TF1 *funcBox = new TF1("funcBox", tbutils::fitFuncBox, xmin, xmax, 3);
	funcBox->SetParameters(8.0, h_residuals->Integral("width"), halfPitch);
	funcBox->FixParameter(3, halfPitch);
	funcBox->SetLineColor(kBlue);
	h_residuals->Fit(funcBox, "+ Q");

	TF1 *funcGaus = new TF1("funcGaus", "gaus", xmin, xmax);
	funcGaus->SetParameters(8.0, h_residuals->Integral("width"), halfPitch);
	funcGaus->FixParameter(3, halfPitch);
	funcGaus->SetLineColor(kGreen);
	h_residuals->Fit(funcGaus, "+ 0 Q");

	makePlot(h_residuals, newHistoName);

	delete[] fileName;
	delete func;
	delete funcBox;
	delete funcGaus;

	return;

}

void Residuals::fitAllClusResiduals(const TbConfig& config,
		const char* analysisName, const char* histoName, TH1 *h_residuals,
		double halfPitch, double xmin, double xmax) {

	char* fileName = config.buildHistName(analysisName, histoName);
	TCanvas* canvas = new TCanvas(TString::Format("can_%s", fileName),
			TString::Format("can_%s", fileName), 10, 20, 700, 500);
	TF1 *multiG = new TF1("multiG", tbutils::multiG, xmin, xmax, 6);
	char newHistoName[150];
	sprintf(newHistoName, "%s-%s", analysisName, histoName);

	double coreMeanGuess = 0.0;
	double outMeanGuess = 0.0;

	double coreWidthGuess = halfPitch * 2.0 / sqrt(12.);
	double outWidthGuess = 6.0 * (halfPitch * 2.0);

	double coreFracGuess = 0.7;
	double coreFracMin = 0.0;
	double coreFracMax = 1.0;

	double normGuess = h_residuals->GetEntries();

	int coreMeanIndex = 0;
	int outMeanIndex = 1;

	int coreWidthIndex = 2;
	int outWidthIndex = 3;

	int coreFracIndex = 4;

	int normIndex = 5;

	multiG->SetParameters(coreMeanGuess, outMeanGuess, coreWidthGuess,
			outWidthGuess, coreFracGuess, normGuess);
	multiG->FixParameter(outMeanIndex, outMeanGuess);
	multiG->SetParLimits(coreFracIndex, coreFracMin, coreFracMax);

	multiG->SetLineColor(kRed);
	h_residuals->Fit(multiG, "Q +");
	//canvas->SaveAs(Form("%s/%s", config.outPath, newHistoName));
	canvas->Write(newHistoName);

	TF1 *func = h_residuals->GetFunction("multiG");
	if (func != NULL) {
		TBALOG(kINFO) << histoName << ": core sigma: "
				<< func->GetParameter(coreWidthIndex) << " outlier sigma: "
				<< func->GetParameter(outWidthIndex) << " core fraction: "
				<< func->GetParameter(coreFracIndex) << " outlier fraction: "
				<< (1. - func->GetParameter(coreFracIndex)) << endl;
	} else {
		TBALOG(kINFO) << histoName << " : fit probably failed" << endl;
	}

	delete[] fileName;
	delete canvas;
	delete multiG;
	delete func;

	return;

}

void Residuals::makePlot(TH1 *h_residuals, char * nameString) {

	TCanvas* canvas = new TCanvas();

	TPad *pad1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0);
	pad1->Draw();
	pad1->cd();
	h_residuals->SetStats(kFALSE);
	h_residuals->Draw();

	TPad *pad2 = new TPad("", "", 0.7, 0.7, 1.0, 1.0);
	pad2->Draw();
	pad2->cd();
	TPaveStats *paramBox = new TPaveStats(0.0, 0.0, 1.0, 1.0);
	paramBox->AddText("Box fit (w/ slope)");
	TF1 *func = h_residuals->GetFunction("func");
	char *nameParam[5] = { (char*) "sigma", (char*) "area", (char*) "x0",
			(char*) "width/2", (char*) "angCoeff" };
	char * line = new char[300];
	for (int i = 0; i < 5; i++) {
		if (func != NULL) {
			sprintf(line, "%s=%10.2e #pm %10.2e", nameParam[i],
					func->GetParameter(i), func->GetParError(i));
		} else {
			sprintf(line, "fit failed?");
		};
		paramBox->AddText(line);
	}
	paramBox->SetFillColor(0);
	paramBox->SetTextColor(kRed);
	paramBox->Draw();

	TPaveStats *paramBox2 = new TPaveStats(0.0, 0.0, 1.0, 1.0);
	paramBox2->AddText("Box fit (w/o slope)");
	TF1 *funcBox = h_residuals->GetFunction("funcBox");
	char *nameParamBox2[3] = { (char*) "sigma", (char*) "area",
			(char*) "width/2" };
	for (int i = 0; i < 3; i++) {
		if (func != NULL) {
			sprintf(line, "%s=%10.2e #pm %10.2e", nameParamBox2[i],
					funcBox->GetParameter(i), funcBox->GetParError(i));
		} else {
			sprintf(line, "fit failed?");
		}
		paramBox2->AddText(line);
	}
	paramBox2->SetFillColor(0);
	paramBox2->SetTextColor(kBlue);
	pad1->cd();
	TPad *pad3 = new TPad("", "", 0.7, 0.5, 1.0, 0.7);
	pad3->Draw();
	pad3->cd();
	paramBox2->Draw();

	TPaveStats *paramBox4 = new TPaveStats(0.0, 0.0, 1.0, 1.0);
	paramBox4->AddText("RMS");
	if (func != NULL) {
		sprintf(line, "%10.2e #pm %10.2e", h_residuals->GetRMS(),
				h_residuals->GetRMSError());
	} else {
		sprintf(line, "fit failed?");
	};
	paramBox4->AddText(line);
	paramBox4->SetFillColor(0);
	paramBox4->SetTextColor(kGreen);
	pad1->cd();
	TPad *pad4 = new TPad("", "", 0.7, 0.15, 1.0, 0.3);
	pad4->Draw();
	pad4->cd();
	paramBox4->Draw();

	//try a gaus if they are close to each other i.e. width id negative
	if (func != NULL) {
		if (func->GetParameter(3) < 4.0) {
			pad1->cd();
			TPaveStats *paramBox3 = new TPaveStats(0.0, 0.0, 1.0, 1.0);
			paramBox3->AddText("Gaus fit");
			TF1 *funcgaus = h_residuals->GetFunction("funcgaus");
			Double_t funcgauspar[3] = { -1, -1, -1 };
			Double_t funcgausparerr[3] = { -1, -1, -1 };
			if (funcgaus) {
				funcgaus->Draw("same");
				for (int i = 0; i < 3; i++) {
					funcgauspar[i] = funcgaus->GetParameter(i);
					funcgausparerr[i] = funcgaus->GetParError(i);
				}

			} else {
				cout << "Problem with fitting Residuals" << endl;
			};
			pad3->cd();
			char *nameParamBox3[3] = { (char*) "area", (char*) "mean",
					(char*) "sigma" };
			for (int i = 0; i < 3; i++) {
				sprintf(line, "%s=%10.2e #pm %10.2e", nameParamBox3[i],
						funcgauspar[i], funcgausparerr[i]);
				sprintf(line, "test");
				paramBox3->AddText(line);
			}
			paramBox3->SetFillColor(0);
			paramBox3->SetTextColor(kGreen);
			pad1->cd();
			TPad *pad4 = new TPad("", "", 0.7, 0.3, 1.0, 0.5);
			pad4->Draw();
			pad4->cd();
			paramBox3->Draw();
		} else {
			cout << "Problem with fitting Residuals" << endl;
		}
	}
	canvas->Write(nameString);
	//canvas->SaveAs(Form("%s/%s", config.outPath, newHistoName));
	delete pad4;
	delete paramBox4;
	delete pad3;
	delete funcBox;
	delete paramBox2;
	delete[] line;
	delete func;
	delete paramBox;
	delete pad2;
	delete pad1;
	delete canvas;
}
