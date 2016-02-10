#include "maxcellresiduals.h"

void MaxCellResiduals::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();
	// hack, but gives reasonable range for both FE-I3 and FE-I4
	int maxToT = (int) floor(2304 / nCols);

	h_resX = new TH1D("", ";X Residuals [#mum]", 800, -1 * pitchX, pitchX);
	h_resY = new TH1D("", ";Y Residuals [#mum]", 800, -1 * pitchX, pitchX);
	h_resXU = new TH1D("", ";Un-cut X Residuals [#mum]", 800, -9000, 9000);
	h_resYU = new TH1D("", ";Un-cut Y Residuals [#mum]", 800, -9000, 9000);
	h_hits = new TH2D("", "", 1000, 0, 60000, 100, -5000, 5000);
	h_hitX = new TH1D("", ";Column; Number of Hits", nCols, 0, nCols);
	h_hitY = new TH1D("", ";Row; Number of Hits", nRows, 0, nRows);
	h_matches = new TH1D("", "", 100, 0, 100000);
	h_sumTot = new TH1D("", ";Max Cell ToT", maxToT, 0, maxToT);
	h_sensorTot = new TH1D("", ";Sensor ToT", maxToT, 0, maxToT);
	h_lv1 = new TH1D("", ";lvl1", 16, 0, 16);
	h_chi2 = new TH1D("", ";Track #chi^{2}", 100, 0, 150);
	m_tracks = 0;
	m_matches = 0;

}

void MaxCellResiduals::event(const TbConfig &config, const Event &event) {

	double pitchX = event.dut->getPitchX();
	double pitchY = event.dut->getPitchY();
	int maxCell = cluster::getMaxTotCell(event.hits);

	for (int ii = 0; ii < event.hits.size(); ii++) {
		h_hitX->Fill(event.hits.at(ii)->col);
		h_hitY->Fill(event.hits.at(ii)->row);
	}
	if (maxCell >= 0) {
		h_lv1->Fill(event.hits.at(0)->lv1);
		h_resYU->Fill(event.hits.at(maxCell)->row * pitchY - event.trackY);
		h_resXU->Fill(event.hits.at(maxCell)->col * pitchX - event.trackX);
		h_hits->Fill(event.hits.at(0)->trig,
				event.hits.at(maxCell)->row * pitchY - event.trackY);
	}
	int maxClu = cluster::getMaxTotCluster(event.clusters);
	if (maxClu >= 0) {
		h_sumTot->Fill(cluster::getSumTot(event.clusters.at(maxClu)));
		h_sensorTot->Fill(cluster::getSumTot(event.hits));
	}

	if (event.fTrack != event::kGood) {
		return;
	}
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	//if( event.ndof <= 4) { return; }

	h_chi2->Fill(event.chi2);
	m_tracks++;
	if (maxCell < 0) {
		return;
	}

	int *a;

	if (event.ndof >= 4) {
		h_resY->Fill(event.hits.at(maxCell)->row * pitchY - event.trackY);
		h_resX->Fill(event.hits.at(maxCell)->col * pitchX - event.trackX);
	}
	h_matches->Fill(event.hits.at(0)->trig);

	for (int ii = 0; ii < event.hits.size(); ii++) {
		if (fabs(event.hits.at(ii)->col * event.dut->getPitchX() - event.trackX)
				> 2.0 * event.dut->getPitchX()) {
			continue;
		}
		if (fabs(event.hits.at(ii)->row * event.dut->getPitchY() - event.trackY)
				> 2.0 * event.dut->getPitchY()) {
			continue;
		}
		m_matches++;
		break;
	}

}

void MaxCellResiduals::finalize(const TbConfig &config) {

	config.drawAndSave(this->name, "lv1", h_lv1);
	config.drawAndSave(this->name, "sumTot", h_sumTot);
	config.drawAndSave(this->name, "sensorTot", h_sensorTot);
	config.drawAndSave(this->name, "resX", h_resX);
	config.drawAndSave(this->name, "resY", h_resY);
	config.drawAndSave(this->name, "row", h_hitX);
	config.drawAndSave(this->name, "col", h_hitY);
	config.drawAndSave(this->name, "resXuncut", h_resXU);
	config.drawAndSave(this->name, "resYuncut", h_resYU);
	config.drawAndSave(this->name, "resVtrig", h_hits);
	config.drawAndSave(this->name, "matchtrig", h_matches);
	config.drawAndSave(this->name, "chi2", h_chi2);
	TBALOG( kINFO ) << "Estimated efficiency : "
			<< double(m_matches) / double(m_tracks) << endl;

}

