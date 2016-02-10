#ifndef TOTCALIB_H
#define TOTCALIB_H

//standard header files
#include <string>
#include <vector>

/*! \brief Base class for all ToT to charge conversion classes.
 *
 * Those not only are responsible for reading
 * the input files whether they are root files or text files, but also for converting the information if necessary,
 * storing it and and proving the methods for making it accessible in a generalized way. The most important
 * method is probably "double q(const int& tot, const int& col = -1, const int&row = -1)".
 */
class ToTCalib {
protected:
	std::string calibration_type;

public:

	std::vector<std::string> calibration_variants;

	// constructor
	ToTCalib();
	// destructor
	~ToTCalib();

	//
	virtual bool addToTCalib(const std::string& filename, const int& iden,
			bool includePerPixel = false) = 0;
	//
	virtual double q(const int& tot, const int& col = -1,
			const int&row = -1) = 0;
	//
	virtual bool hasPerPixelCalib() = 0;
	//
	virtual bool hasCalib() = 0;
	//
	virtual std::string getCalibrationType(){
		return calibration_type;
	};

};


#endif //TOTCALIB_H
