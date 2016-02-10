#include "beamprofile.h"

void BeamProfile::init(TbConfig &config) {

	const DUT* dut = config.getDut(this->iden);
	double pitchX = dut->getPitchX();
	double pitchY = dut->getPitchY();
	int nCols = dut->getNcols();
	int nRows = dut->getNrows();

	h_beamProfile = new TH2D("", ";Track X [#mum];Track Y [#mum]", 200, -15000,
			15000, 200, -15000, 15000);
	h_beamProfileZoom = new TH2D("", ";Track X [#mum];Track Y [#mum]", 50,
			-1 * pitchX, nCols * pitchX, 50, -1 * pitchY, nRows * pitchY);
	h_beamProfileMasked = new TH2D("", ";Track X [#mum];Track Y [#mum]", nCols,
			0, nCols * pitchX, nRows, 0, nRows * pitchY);
	h_dutHit = new TH2D("", ";Track X [#mum];Track Y [#mum]", 80, -5000, 15000,
			80, -5000, 15000);
	h_chi2 = new TH1D("", ";Track #chi^{2}", 100, 0, 300);
	h_chi2Zoom = new TH1D("", ";Track #chi^{2}", 100, 0, 30);

	h_angle = new TH2D("", ";dx/dz;dy/dz", 400, -0.004, 0.004, 400, -0.004,
			0.004);
	h_angle_normalcuts = new TH2D("", ";dx/dz;dy/dz", 400, -0.004, 0.004, 400,
			-0.004, 0.004);
	h_beamProfile_noAngleCut = (TH2D*) h_beamProfile->Clone(""); //Same binning etc.

	h_corr_PosVsAngle_XX = new TH2D("", "", 400, -15000, 15000, 400, -0.004,
			0.004);
	h_corr_PosVsAngle_XY = (TH2D*) h_corr_PosVsAngle_XX->Clone("");
	h_corr_PosVsAngle_YX = (TH2D*) h_corr_PosVsAngle_XX->Clone("");
	h_corr_PosVsAngle_YY = (TH2D*) h_corr_PosVsAngle_XX->Clone("");
	h_corr_PosVsAngle_XX->SetTitle(";x [#mum];dx/dz");
	h_corr_PosVsAngle_XY->SetTitle(";x [#mum];dy/dz");
	h_corr_PosVsAngle_YX->SetTitle(";y [#mum];dx/dz");
	h_corr_PosVsAngle_YY->SetTitle(";y [#mum];dy/dz");

}

void BeamProfile::event(const TbConfig &config, const Event &event) {

	TBALOG(kDEBUG2) << "Top of BeamProfile::event" << endl;

	h_chi2->Fill(event.chi2);
	h_chi2Zoom->Fill(event.chi2);

	//DON'T check angleCuts, only Chi2
	if (event.fBase == event::kGood and event.fTrackChi2 == event::kGood) {
		h_angle->Fill(event.dxdz, event.dydz);
		//cout << "dxdz=" << event.dxdz << " dydz=" << event.dydz << endl; //Useful for finding the scale
		h_beamProfile_noAngleCut->Fill(event.trackX, event.trackY);
		h_corr_PosVsAngle_XX->Fill(event.trackX, event.dxdz);
		h_corr_PosVsAngle_XY->Fill(event.trackX, event.dydz);
		h_corr_PosVsAngle_YX->Fill(event.trackY, event.dxdz);
		h_corr_PosVsAngle_YY->Fill(event.trackY, event.dydz);
	}

	if (event.fTrack != event::kGood) {
		return;
	}

	h_beamProfile->Fill(event.trackX, event.trackY);
	h_angle_normalcuts->Fill(event.dxdz, event.dydz);

	if (event.hits.size() > 0) {
		h_dutHit->Fill(event.trackX, event.trackY);
	}

	if (event.fTrackRegion == event::kGood) {
		h_beamProfileMasked->Fill(event.trackX, event.trackY);
	}

	if (event.trackX > -1.0 * event.dut->getPitchX()
			&& event.trackX < event.dut->getNcols() * event.dut->getPitchX()) {
		if (event.trackY > -1.0 * event.dut->getPitchY()
				&& event.trackY
						< event.dut->getNrows() * event.dut->getPitchY()) {
			h_beamProfileZoom->Fill(event.trackX, event.trackY);
		}
	}

}

void BeamProfile::finalize(const TbConfig &config) {

	config.drawAndSave(name, "beamprofile", "cont4z", h_beamProfile);
	config.drawAndSave(name, "beamprofilezoom", "colz", h_beamProfileZoom);
	config.drawAndSave(name, "beamprofilemasked", "colz", h_beamProfileMasked);
	config.drawAndSave(name, "beamprofile_noAngleCut", "cont4z",
			h_beamProfile_noAngleCut);
	config.drawAndSave(name, "duthit", "cont4z", h_dutHit);
	config.drawAndSave(name, "chi2", h_chi2);
	config.drawAndSave(name, "chi2Zoom", h_chi2Zoom);
	config.drawAndSave(name, "angle", "colz", h_angle);
	config.drawAndSave(name, "angle_normalcuts", "colz", h_angle_normalcuts);
	config.drawAndSave(name, "corr_PosVsAngle_XX", h_corr_PosVsAngle_XX);
	config.drawAndSave(name, "corr_PosVsAngle_XY", h_corr_PosVsAngle_XY);
	config.drawAndSave(name, "corr_PosVsAngle_YX", h_corr_PosVsAngle_YX);
	config.drawAndSave(name, "corr_PosVsAngle_YY", h_corr_PosVsAngle_YY);

	TH1D* h_trackX = (TH1D*) h_beamProfile->ProjectionX(" ");
	config.drawAndSave(name, "trackX", h_trackX);
	TH1D* h_trackY = (TH1D*) h_beamProfile->ProjectionY(" ");
	config.drawAndSave(name, "trackY", h_trackY);

	TH1D* h_angle_X = (TH1D*) h_angle->ProjectionX(" ");
	config.drawAndSave(name, "angle_dxdz", h_angle_X);
	TH1D* h_angle_Y = (TH1D*) h_angle->ProjectionY(" ");
	config.drawAndSave(name, "angle_dydz", h_angle_Y);

	TH1D* h_trackX_noAngleCut = (TH1D*) h_beamProfile_noAngleCut->ProjectionX(
			" ");
	config.drawAndSave(name, "trackX_noAngleCut", h_trackX_noAngleCut);
	TH1D* h_trackY_noAngleCut = (TH1D*) h_beamProfile_noAngleCut->ProjectionY(
			" ");
	config.drawAndSave(name, "trackY_noAngleCut", h_trackY_noAngleCut);

}
