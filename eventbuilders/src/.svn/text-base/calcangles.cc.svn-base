#include "calcangles.h"
#include "dut.h"

/// Calculates tilt of the sensor 
/// Documentation on algorithm can be found in tbmon/doc
/// Well, at least soon. TODO: Make documentation
void CalcAngles::initRun(TbConfig &config) {
	double pi = atan2(0.0f, -1.0f);
	for (list<DUT*>::iterator it = config.dutList.begin();
			it != config.dutList.end(); it++) {
		DUT *dut = *it;
		dut->anglesCalculated = false;

		if (dut->hasAnglesFromReco) {
			double cosPhiGear = cos(dut->anglePhiFromGEAR);
			double sinPhiGear = sin(dut->anglePhiFromGEAR);
			double cosEtaGear = cos(dut->angleEtaFromGEAR);
			double sinEtaGear = sin(dut->angleEtaFromGEAR);
			double cosPhiAlign = cos(dut->anglePhiFromAlignment);
			double sinPhiAlign = sin(dut->anglePhiFromAlignment);
			double cosEtaAlign = cos(dut->anglePhiFromAlignment);
			double sinEtaAlign = sin(dut->anglePhiFromAlignment);

			double denom = cosEtaAlign
					* (cosPhiGear * cosEtaGear * cosPhiAlign
							- sinPhiGear * sinPhiAlign)
					+ (cosPhiGear * sinEtaGear * sinEtaAlign);

			TBALOG(kDEBUG2) << "DUT " << dut->getDUTid() << " denom: " << denom
					<< endl;
			if (denom == 0) {
				TBALOG(kERROR) << "DUT " << dut->getDUTid()
						<< ":Could not compute angles. Division by zero."
						<< endl;
				continue;
			}

			double phiNumerator = cosEtaAlign
					* (cosPhiGear * sinPhiAlign
							+ sinPhiGear * cosEtaGear * cosPhiAlign)
					- sinPhiGear * sinEtaGear * sinEtaAlign;
			TBALOG(kDEBUG2) << "DUT " << dut->getDUTid() << ": phiNumerator: "
					<< phiNumerator << endl;

			double etaNumerator = -sinEtaGear * cosEtaAlign * cosPhiAlign
					- cosEtaGear * sinEtaAlign;
			TBALOG(kDEBUG2) << "DUT " << dut->getDUTid() << ": etaNumerator: "
					<< etaNumerator << endl;

			dut->anglePhiCalculated = atan(phiNumerator / denom);
			dut->angleEtaCalculated = atan(etaNumerator / denom);
			dut->anglesCalculated = true;
			TBALOG(kINFO) << "DUT " << dut->getDUTid()
					<< ": Calculated angle Phi: "
					<< dut->anglePhiCalculated * 180 / pi << " deg" << endl;
			TBALOG(kINFO) << "DUT " << dut->getDUTid()
					<< ": Calculated angle Eta: "
					<< dut->angleEtaCalculated * 180 / pi << " deg" << endl;
		} else {
			TBALOG(kERROR) << "DUT " << dut->getDUTid()
					<< ": No angles found in DUT object" << endl;
		}
	}
}
