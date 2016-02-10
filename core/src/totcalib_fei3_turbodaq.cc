#include "totcalib_fei3_turbodaq.h"

const int ToTCalib_FEI3_TurboDAQ::rows = 160;
const int ToTCalib_FEI3_TurboDAQ::cols = 18;

double ToTCalib_FEI3_TurboDAQ::q(const int& tot, const int& col,const int&row) {
	if (calibpp[index(col,row)]==NULL) {
		std::cerr << "Error with charge conversion, missing calibration values for col " << col << " row " << row << ". Returning zero" << std::endl;
		return -1;
	}
	return calibpp[index(col,row)]->calB / (tot - calibpp[index(col,row)]->calA) - calibpp[index(col,row)]->calC;
}

bool ToTCalib_FEI3_TurboDAQ::addToTCalib(const std::string& filename, const int& iden, bool includePerPixel) {
	int chip,column,row;
	float scanvar,parA,parB,parC;
	char line[255];
	int success;
	fstream calibfile;
	bool readin_successfull;

	readin_successfull = true;

	calibfile.open(filename.c_str());
	if (! calibfile.is_open()) {
		std::cerr << "Calibfile is not open" << std::endl;
		readin_successfull = false;
	}
	else while (!calibfile.eof()) {
		calibfile.getline(line,255);
        	success = sscanf(line,"%d %d %d %f %f %f %f",&chip,&column,&row,&scanvar,&parA,&parB,&parC);
		while (calibfile.peek() < 33 && !calibfile.eof()) {calibfile.ignore(1);} //throw away any stray end of line characters
		if (chip+10!=iden) continue; //check have correct chip ID
		if (success==7) {calibpp[index(column,row)] =  new calib(parA,parB,parC);}
		else {
			std::cerr << "Error with read in of calibratory values for charge conversion. Results may be somewhat dubious" << std::endl;
			readin_successfull = false;
		}
	}
	return readin_successfull;
}

