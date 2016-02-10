#include "totcalib.h"


ToTCalib::ToTCalib() {
	calibration_variants.push_back("fei3_usbpix_converted");
	calibration_variants.push_back("fei3_turbodaq_textfile");
	calibration_variants.push_back("fei3_turbodaq_converted_from_tot_file");
	calibration_variants.push_back("undefined");
	calibration_type = "undefined";
};

ToTCalib::~ToTCalib() {}

