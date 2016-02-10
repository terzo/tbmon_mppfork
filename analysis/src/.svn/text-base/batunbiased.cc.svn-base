#include "batunbiased.h"

using namespace std;

const int BatUnbiased::batIdens[3] = { 81, 83, 86 };

void BatUnbiased::init(TbConfig& config) {
	if (config.logLevel == kDEBUG3) {
		cout << "[ batunbiased ]; In init()" << endl;
	}

	h_residuals = new TH1D**[numBats];
	h_residuals_zoom = new TH1D**[numBats];

	h_etaCorr = new TH1D**[numBats];
	h2_residualsVsEtaCorr = new TH2D**[numBats];

	h_trackVsTruth = new TH1D**[numBats];
	h_hitVsTruth = new TH1D**[numBats];

	h2_trueresidualVSinterstrip_truehit = new TH2D**[numBats];
	h2_trueresidualVsEta = new TH2D**[numBats];
	h2_trueresidualVsEtaCorr = new TH2D**[numBats];

	h2_stripVstruth = new TH2D***[numBats];

	for (int i = 0; i < numBats; i++) {
		h_residuals[i] = new TH1D*[2];
		h_residuals_zoom[i] = new TH1D*[2];

		h_etaCorr[i] = new TH1D*[2];
		h2_residualsVsEtaCorr[i] = new TH2D*[2];

		h_trackVsTruth[i] = new TH1D*[2];
		h_hitVsTruth[i] = new TH1D*[2];

		h2_trueresidualVSinterstrip_truehit[i] = new TH2D*[2];
		h2_trueresidualVsEta[i] = new TH2D*[2];
		h2_trueresidualVsEtaCorr[i] = new TH2D*[2];

		h2_stripVstruth[i] = new TH2D**[2];

		for (int j = 0; j <= 1; j++) {
			char title[100];
			char name[100];

			/*
			 * 2nd index is BAT LAYER
			 *
			 * Old naming convention, should be removed:
			 *  "Y" = straight up                                       = layer0
			 *  "X" = horizontal, pointing left when looking downstream = layer1
			 */

			sprintf(title, "BAT-%i, Unbiased residuals (layer %i)", batIdens[i],
					j);
			sprintf(name, "BAT-%i_residuals_%i", batIdens[i], j);
			h_residuals[i][j] = new TH1D(name, title, 400, -2000, 2000);
			h_residuals[i][j]->SetXTitle("Track - hit [#mu m]");

			sprintf(title, "BAT-%i, Unbiased residuals (zoom) (layer %i)",
					batIdens[i], j);
			sprintf(name, "BAT-%i_residuals_zoom_%i", batIdens[i], j);
			h_residuals_zoom[i][j] = new TH1D(name, title, 400, -50, 50);
			h_residuals_zoom[i][j]->SetXTitle("(track - hit) [#mu m]");

			sprintf(title,
					"BAT-%i, Unbiased residuals vs. eta-corrected hit position estimate (layer %i)",
					batIdens[i], j);
			sprintf(name, "BAT-%i_residualsVsEtaCorr_%i", batIdens[i], j);
			h2_residualsVsEtaCorr[i][j] = new TH2D(name, title, 100, 0.0, 1.0,
					400, -50, 50);
			h2_residualsVsEtaCorr[i][j]->SetYTitle("(track - hit) [#mu m]");
			h2_residualsVsEtaCorr[i][j]->SetXTitle("Eta-correction");

			sprintf(title,
					"BAT-%i, Eta-corrected interstrip hit estimator (layer %i)",
					batIdens[i], j);
			sprintf(name, "BAT-%i_etaCorr_%i", batIdens[i], j);
			h_etaCorr[i][j] = new TH1D(name, title, 100, 0.0, 1.0);
			h_etaCorr[i][j]->SetXTitle("Eta-correction");

			if (config.isSimulation) {
				sprintf(title, "BAT-%i, Track Vs Truth (layer %i)", batIdens[i],
						j);
				sprintf(name, "BAT-%i_TrackVsTruth_%i", batIdens[i], j);
				h_trackVsTruth[i][j] = new TH1D(name, title, 400, -50, 50);
				h_trackVsTruth[i][j]->SetXTitle("(track - truth) [#mu m]");

				sprintf(title, "BAT-%i, Hit Vs Truth (layer %i)", batIdens[i],
						j);
				sprintf(name, "BAT-%i_HitVsTruth_%i", batIdens[i], j);
				h_hitVsTruth[i][j] = new TH1D(name, title, 400, -50, 50);
				h_hitVsTruth[i][j]->SetXTitle("(hit - truth) [#mu m]");

				sprintf(title,
						"BAT-%i, Truth-residual Vs interstrip-truth (layer %i)",
						batIdens[i], j);
				sprintf(name, "BAT-%i_trueresidualVSinterstrip_truehit_%i",
						batIdens[i], j);
				h2_trueresidualVSinterstrip_truehit[i][j] = new TH2D(name,
						title, 400, 0, 1, 400, -50, 50);
				h2_trueresidualVSinterstrip_truehit[i][j]->SetXTitle(
						"True hit position modulo strip pitch");
				h2_trueresidualVSinterstrip_truehit[i][j]->SetYTitle(
						"Truth - hit estimator [#mu m]");

				sprintf(title, "BAT-%i, Truth-residual Vs eta (layer %i)",
						batIdens[i], j);
				sprintf(name, "BAT-%i_trueresidualVSeta_%i", batIdens[i], j);
				h2_trueresidualVsEta[i][j] = new TH2D(name, title, 400, 0.0,
						1.0, 400, -50, 50);
				h2_trueresidualVsEta[i][j]->SetXTitle("3-strip #eta [#mu m]");
				h2_trueresidualVsEta[i][j]->SetYTitle(
						"Truth - hit estimator [#mu m]");

				sprintf(title,
						"BAT-%i, Truth-residual Vs eta-corrected hit position estimate (layer %i)",
						batIdens[i], j);
				sprintf(name, "BAT-%i_trueresidualVSetaCorr_%i", batIdens[i],
						j);
				h2_trueresidualVsEtaCorr[i][j] = new TH2D(name, title, 400, 0.0,
						1.0, 400, -50, 50);
				h2_trueresidualVsEtaCorr[i][j]->SetXTitle("#eta-correction");
				h2_trueresidualVsEtaCorr[i][j]->SetYTitle(
						"Truth - hit estimator [#mu m]");

				//Correlation histos
				h2_stripVstruth[i][j] = new TH2D*[2]; //x,y
				for (int k = 0; k <= 1; k++) {
					sprintf(title,
							"BAT-%i, strip (layer %i) Vs Truth (tbmon %c)",
							batIdens[i], j, k == 0 ? 'x' : 'y');
					sprintf(name, "BAT-%i_strip_vs_truth_(%iv%c)", batIdens[i],
							j, k == 0 ? 'x' : 'y');
					h2_stripVstruth[i][j][k] = new TH2D(name, title, 640, 0,
							640, 100, -32000, 32000);
				}
			}
		}
	}
}

void BatUnbiased::event(const TbConfig &config, const Event &event) {
	//Cuts
	if (event.fTrack == event::kBad || event.fBase == event::kBad)
		return;
	if (event.track == NULL)
		return; //Not applicable

	//Analyze all BAT's. Just set a single iden
	for (int i = 0; i < numBats; i++) {
		int currentIden = batIdens[i];

		//Find the trackParams
		TrackParams* params = NULL;
		for (int iPar = 0; params == NULL && iPar < event.track->nTrackParams;
				iPar++) {
			params = (TrackParams*) event.track->trackParams[iPar];
			if (params->iden == currentIden)
				break; //Found
			params = NULL; //Not found
		}
		if (params == NULL) {
			cout << "[ batUnbiased::event ]; TrackParams for iden"
					<< currentIden << " not found!" << endl;
			continue; //Skip this BAT in this event
		}

		//Track vs. hit
		//Loop over batClusters, Fill when iden/layer matches
		for (int iCluster = 0; iCluster < event.track->nBatCluster;
				iCluster++) {
			BatCluster* cluster =
					(BatCluster*) event.track->batCluster[iCluster];
			if (cluster->iden == currentIden) {
				double residual = params->params[cluster->claye ? 0 : 1] * 1000
						- cluster->ccorr * 50;
				h_residuals[i][cluster->claye]->Fill(residual);
				h_residuals_zoom[i][cluster->claye]->Fill(residual);

				double etaCorr = fmod((double) cluster->ccorr, 1.0);
				h2_residualsVsEtaCorr[i][cluster->claye]->Fill(etaCorr,
						residual);
				h_etaCorr[i][cluster->claye]->Fill(etaCorr);
			}
		}

		//Simulation-specific (need truth info)
		if (config.isSimulation) {
			//Get truth info
			map<int, string>::const_iterator m_it =
					config.simIdenNameMap_all.find(currentIden);
			string simModName =
					(m_it == config.simIdenNameMap_all.end()) ?
							"noDevice" : m_it->second;

			for (vector<simTruthHit*>::iterator truthHit =
					event.simData->allTruthHits->begin();
					truthHit != event.simData->allTruthHits->end();
					truthHit++) {
				if ((*truthHit)->planeID == simModName) {
					//We have a valid truthHit (there might be more than one): Analyze!

					//Transform truthHit's posLocal into coordinate system as used by tbreco
					simThreeVector posLocal = (*truthHit)->posLocalRaw;
					double tmpPos = posLocal.data[1];
					posLocal.data[1] = posLocal.data[0] * 1000
							+ (640 / 2.0 - 0.5) * 50;
					posLocal.data[0] = tmpPos * 1000 + (640 / 2.0 - 0.5) * 50;

					//Track vs. truth
					//Loop over layer (y,x)
					for (int layer = 0; layer <= 1; layer++) {
						h_trackVsTruth[i][layer]->Fill(
								params->params[layer] * 1000
										- posLocal.data[layer ? 0 : 1]);
					}

					//Hit vs. Truth, True residual vs. interstrip true hit position, and correlation histo
					//Loop over batClusters, fill when iden matches
					for (int iCluster = 0; iCluster < event.track->nBatCluster;
							iCluster++) {
						BatCluster* cluster =
								(BatCluster*) event.track->batCluster[iCluster];
						if (cluster->iden == currentIden) {
							double trueRes = posLocal.data[cluster->claye]
									- cluster->ccorr * 50;
							double interTrueHit = fmod(
									posLocal.data[cluster->claye], 50.0) / 50.0;
							h_hitVsTruth[i][cluster->claye]->Fill(trueRes);

							h2_trueresidualVSinterstrip_truehit[i][cluster->claye]->Fill(
									interTrueHit, trueRes);
							double eta = fmod((double) cluster->cstrp, 1.0);
							//cout << currentIden << " " << cluster->claye << " " << cluster->ccorr << " " <<cluster->cstrp << endl;
							h2_trueresidualVsEta[i][cluster->claye]->Fill(eta,
									trueRes);

							double etaCorr = fmod((double) cluster->ccorr, 1.0);
							h2_trueresidualVsEtaCorr[i][cluster->claye]->Fill(
									etaCorr, trueRes);

							for (int p = 0; p <= 1; p++) {
								h2_stripVstruth[i][cluster->claye][p]->Fill(
										cluster->ccorr, posLocal.data[p]);
							}
						}
					}
				}
			} //End loop over truthHits
		} //End simulation analysis
	} //End loop over BATs
}

void BatUnbiased::finalize(const TbConfig &config) {
	for (int i = 0; i < numBats; i++) {
		for (int j = 0; j <= 1; j++) {
			config.drawToFile(this->name, h_residuals[i][j]->GetName(),
					h_residuals[i][j]);
			config.saveToFile(this->name, h_residuals[i][j]->GetName(),
					h_residuals[i][j]);

			config.drawToFile(this->name, h_residuals_zoom[i][j]->GetName(),
					h_residuals_zoom[i][j]);
			config.saveToFile(this->name, h_residuals_zoom[i][j]->GetName(),
					h_residuals_zoom[i][j]);

			config.drawToFile(this->name,
					h2_residualsVsEtaCorr[i][j]->GetName(),
					h2_residualsVsEtaCorr[i][j]);
			config.saveToFile(this->name,
					h2_residualsVsEtaCorr[i][j]->GetName(),
					h2_residualsVsEtaCorr[i][j]);

			config.drawAndSave(this->name, h_etaCorr[i][j]->GetName(),
					h_etaCorr[i][j]);

			if (config.isSimulation) {
				config.drawToFile(this->name, h_trackVsTruth[i][j]->GetName(),
						h_trackVsTruth[i][j]);
				config.saveToFile(this->name, h_trackVsTruth[i][j]->GetName(),
						h_trackVsTruth[i][j]);

				config.drawToFile(this->name, h_hitVsTruth[i][j]->GetName(),
						h_hitVsTruth[i][j]);
				config.saveToFile(this->name, h_hitVsTruth[i][j]->GetName(),
						h_hitVsTruth[i][j]);

				config.saveToFile(this->name,
						h2_trueresidualVSinterstrip_truehit[i][j]->GetName(),
						h2_trueresidualVSinterstrip_truehit[i][j]);
				config.drawToFile(this->name,
						h2_trueresidualVSinterstrip_truehit[i][j]->GetName(),
						h2_trueresidualVSinterstrip_truehit[i][j]);

				config.saveToFile(this->name,
						h2_trueresidualVsEta[i][j]->GetName(),
						h2_trueresidualVsEta[i][j]);
				config.drawToFile(this->name,
						h2_trueresidualVsEta[i][j]->GetName(),
						h2_trueresidualVsEta[i][j]);

				config.drawAndSave(this->name,
						h2_trueresidualVsEtaCorr[i][j]->GetName(),
						h2_trueresidualVsEtaCorr[i][j]);

				for (int k = 0; k <= 1; k++) {
					config.drawToFile(this->name,
							h2_stripVstruth[i][j][k]->GetName(),
							h2_stripVstruth[i][j][k]);
					config.saveToFile(this->name,
							h2_stripVstruth[i][j][k]->GetName(),
							h2_stripVstruth[i][j][k]);
				}
			}
		}
	}
}
