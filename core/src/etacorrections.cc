#include "etacorrections.h"

void EtaCorrections::addEcorr(const char* fileName, int nbins) {

	EtaCorr corr;
	corr.name = (char*) fileName;
	corr._nbins = nbins;

	ifstream stream;
	stream.open(fileName);
	if (not stream.is_open()) {
		cerr << "Unable to open file " << fileName << endl;
		exit(-1);
	}
	stream >> corr.firstRun;
	stream >> corr.lastRun;
	corr.ecorrY.resize(corr._nbins);
	double tmp;
	for (int ii = 0; ii < corr._nbins; ii++) {
		stream >> tmp;
		corr.ecorrY.at(ii) = tmp;
	}
	corr.ecorrX.resize(corr._nbins);
	for (int ii = 0; ii < corr._nbins; ii++) {
		stream >> tmp;
		corr.ecorrX.at(ii) = tmp;
	}
	allCorrs.push_back(corr);

	stream >> tmp;

	if (not stream.eof()) {
		cout << "Malformed eta calib file, all values not read." << endl;
	}
}

void EtaCorrections::initRun(int currentRun) {
	bool gotIt(false);
	for (vector<EtaCorr>::iterator it = allCorrs.begin(); it != allCorrs.end();
			it++) {
		if (currentRun < (*it).firstRun) {
			continue;
		}
		if (currentRun > (*it).lastRun) {
			continue;
		}
		eCorrs = (*it);
		gotIt = true;
		break;
	}
	//cout << "Using calibs from file: " << eCorrs.name << endl;
	if (not gotIt) {
		cout << "Eta corrections not found for run " << currentRun << endl;
	}
}

double EtaCorrections::getX(double chargeCorrected) {
	int bin = (int) floor(chargeCorrected * eCorrs._nbins);
	return bin < eCorrs.ecorrX.size() ? eCorrs.ecorrX.at(bin) : 0.0;
}

double EtaCorrections::getY(double chargeCorrected) {
	int bin = (int) floor(chargeCorrected * eCorrs._nbins);
	//cout << "Ecorr bin: " << bin << "Val: " << eCorrs.ecorrY.at(bin) <<endl;
	return bin < eCorrs.ecorrY.size() ? eCorrs.ecorrY.at(bin) : 0.0;
}
