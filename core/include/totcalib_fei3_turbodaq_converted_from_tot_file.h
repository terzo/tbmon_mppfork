#ifndef ROOT_TOTCALIB_H
#define ROOT_TOTCALIB_H

// standard header files
#include <assert.h>
#include <map>
#include <stdlib.h>

// root header files
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"

// tbmon header files
#include "totcalib.h"

/*! \brief Class reads calibration for FE-I3 from root files prepared by script in the directory totcalib from TurboDAQ *.tot files.
 *
 *
 */
class ToTCalib_FEI3_TurboDaq_converted_from_ToT_file: public ToTCalib {

public:
	ToTCalib_FEI3_TurboDaq_converted_from_ToT_file();
	~ToTCalib_FEI3_TurboDaq_converted_from_ToT_file();

	bool addToTCalib(const std::string& filename, const int& iden,
			bool includePerPixel = false);
	bool hasPerPixelCalib() {
		return calibpp.empty() ? false : true;
	}
	bool hasCalib() {
		return calib == NULL ? false : true;
	}
	double q(const int& tot, const int& col = -1, const int&row = -1) {
		double electrons =
				(col < 0 || row < 0) ? charge(tot) : chargepp(tot, col, row);
		electrons *= chargeCorrection();
		return electrons;
	}
	void chargeCorrection(const double& c) {
		correction = c;
	}
	double chargeCorrection() {
		return correction;
	}
	void checkRange(const int& tot);

private:
	TF1* calib;
	std::map<int, TF1*> calibpp;
	static const int rows;
	static const int cols;

	int index(const int& col, const int& row) {
		return rows * row + col;
	}
	double chargepp(const int& tot, const int& col, const int&row);
	double charge(const int& tot);
	double correction;
};

#endif //ROOT_TOTCALIB_H
