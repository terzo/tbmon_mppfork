#include "checktrack.h"

void CheckTrack::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	h_chi2 = new TH1D("", ";Track #chi^{2}", 100, 0, 100);
	h_badRegion = new TH2D("", ";Column;Row", nCols, 0, nCols, nRows, 0, nRows);

}

void CheckTrack::event(const TbConfig &config, const Event &event) {

	h_chi2->Fill(event.chi2);
	if (event.fTrackRegion != event::kGood) {
		if (event.fTrackCentralRegion == event::kGood) {
			TBALOG(kDEBUG2) << "Masked at "
					<< tbutils::getCol(event.trackX, event) << ", "
					<< tbutils::getRow(event.trackY, event) << endl;
			h_badRegion->Fill(tbutils::getCol(event.trackX, event),
					tbutils::getRow(event.trackY, event));
		}
	}

}

void CheckTrack::finalize(const TbConfig &config) {

	config.drawAndSave(this->name, "chi2", h_chi2);
	config.drawAndSave(this->name, "badRegion", "col", h_badRegion);

}
