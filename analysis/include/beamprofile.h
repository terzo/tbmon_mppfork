#ifndef BEAMPROFILE_H
#define BEAMPROFILE_H

// root header files
#include <TH2D.h>

// tbmon header files
#include "event.h"
#include "tbanalysis.h"

/*! \brief shows where and how the beam hit the sensors; track maps, track chi2, beam angle
distributions
 *
 */
class BeamProfile: public TbAnalysis {
private:
	TH2D* h_beamProfile;
	TH2D* h_beamProfileZoom;
	TH2D* h_beamProfileMasked;
	TH1D* h_chi2;
	TH1D* h_chi2Zoom;
	TH2D* h_dutHit;
	TH2D* h_angle;
	TH2D* h_angle_normalcuts;
	TH2D* h_beamProfile_noAngleCut;

	TH2D* h_corr_PosVsAngle_XX;
	TH2D* h_corr_PosVsAngle_XY;
	TH2D* h_corr_PosVsAngle_YX;
	TH2D* h_corr_PosVsAngle_YY;

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //BEAMPROFILE_H
