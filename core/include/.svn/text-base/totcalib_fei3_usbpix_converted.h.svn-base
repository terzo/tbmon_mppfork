#ifndef TOTCALIB_FEI3_USBPIX_CONVERTED_H
#define TOTCALIB_FEI3_USBPIX_CONVERTED_H

// standard header files
#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>

// root header files
#include "TDirectory.h"
#include "TFile.h"

// tbmon header files
#include "totcalib.h"

/*! \brief Class reads calibration for FE-I3 from text files converted by the (e.g.) the ABCtotal.py script from USBpix created root files to text files. Those files are similar to native TurboDAQ TOTcalib scans but the conversion equation has changed slightly.
 *
 *
 */
class ToTCalib_FEI3_USBpix_converted: public ToTCalib {

private:

	/*! \brief contains calibration constants for one pixel in class ToTCalib_FEI3_USBpix_converted
	 *
	 */
	struct calib {
		float calA; ///< calibration constant A
		float calB; ///< calibration constant B
		float calC; ///< calibration constant C
		// sets calibration constants for one pixel in class ToTCalib_FEI3_USBpix_converted
		calib(float parA, float parB, float parC) {
			calA = parA;
			calB = parB;
			calC = parC;
		}
	};

	static const int rows; ///< number of rows - fixed and not read from DUT as this class only deals with FE-I3
	static const int cols; ///< number of columns - fixed and not read from DUT as this class only deals with FE-I3

	// serialization method for pixels in FE-I3
	int index(const int& col, const int& row) {
		return rows * row + col;
	}

	// calibration constants for one DUT/FE, serialized by index(col, row) method
	std::map<int, calib*> calibpp;

public:

	ToTCalib_FEI3_USBpix_converted() {
		calibration_type = "fei3_usbpix_converted";
	}
	;

	~ToTCalib_FEI3_USBpix_converted() {
	}
	;

	/* add calibration for one DUT
	 * returns true if successfull, otherwise false
	 */
	bool addToTCalib(const std::string& filename, const int& iden,
			bool includePerPixel = false);
	/* turns tot into charge using calibration constants of a specific pixel,
	 * returns -1 if calibration for this pixel does not exist and throws error message on
	 * std:cerr
	 */
	double q(const int& tot, const int& col = -1, const int&row = -1);
	/* basically returns whether calibrations has been (successfully) already loaded for this DUT as well,
	 * as in this case it has always a pixel calibration
	 */
	bool hasPerPixelCalib() {
		return calibpp.empty() ? false : true;
	}
	// returns whether calibrations has been (successfully) already loaded for this DUT
	bool hasCalib() {
		return calibpp.empty() ? false : true;
	}
};

#endif //TOTCALIB_FEI3_USBPIX_CONVERTED
