#include "correlations.h"

void Correlation::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	h_colvx = new TH2D("", ";Raw Hits Position [col];Track X [#mum]", nCols, 0,
			nCols, 64, -9000, 9000);
	h_colvy = new TH2D("", ";Raw Hits Position [col];Track Y [#mum]", nCols, 0,
			nCols, 64, -9000, 9000);
	h_rowvx = new TH2D("", ";Raw Hits Position [row];Track X [#mum]", nRows, 0,
			nRows, 64, -9000, 9000);
	h_rowvy = new TH2D("", ";Raw Hits Position [row];Track Y [#mum]", nRows, 0,
			nRows, 64, -9000, 9000);
	h_occup = new TH2D("", ";Column;Row", nCols, 0, nCols, nRows, 0, nRows);

}

void Correlation::event(const TbConfig &config, const Event &event) {

	double trkX(0), trkY(0);
	trkX = event.trackX;
	trkY = event.trackY;
	TBALOG(kDEBUG3) << "trkX = " << trkX << ", trkY = " << trkY << endl;

	for (int ii = 0; ii < event.rawHits.size(); ii++) {
		h_colvx->Fill(event.rawHits.at(ii)->col, trkX);
		h_colvy->Fill(event.rawHits.at(ii)->col, trkY);
		h_rowvx->Fill(event.rawHits.at(ii)->row, trkX);
		h_rowvy->Fill(event.rawHits.at(ii)->row, trkY);
		h_occup->Fill(event.rawHits.at(ii)->col, event.rawHits.at(ii)->row);
	}

}

void Correlation::finalize(const TbConfig &config) {

	config.drawAndSave(this->name, "colvx", h_colvx);
	config.drawAndSave(this->name, "colvy", h_colvy);
	config.drawAndSave(this->name, "rowvx", h_rowvx);
	config.drawAndSave(this->name, "rowvy", h_rowvy);
	config.drawAndSave(this->name, "occup", "colz", h_occup);

}
