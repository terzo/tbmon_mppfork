#ifndef ETACORRECTIONS_H
#define ETACORRECTIONS_H

// standard header files
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

using namespace std;
/*! \brief Handles read-in from file and storage of eta correction values. Only instantiated in DUT.
 *
 */
class EtaCorrections {
public:
	/*! \brief Contains eta correction data for one DUT.
	 *
	 */
	struct EtaCorr {
		char* name;
		int firstRun;
		int lastRun;
		vector<double> ecorrX;
		vector<double> ecorrY;
		int _nbins;
	};
private:
	vector<EtaCorr> allCorrs;
public:
	EtaCorr eCorrs;
	void initRun(int currentRun);
	void addEcorr(const char* fileName, int nbins = 100);
	double getX(double chargeCorrected);
	double getY(double chargeCorrected);
};

#endif //ETACORRECTIONS_H
