#include "totcalib_fei3_turbodaq_converted_from_tot_file.h"

using namespace std;

const int ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::rows = 160;
const int ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::cols = 18;

double conv(double* x, double* par) {
	//From Sara S.
	//Fit to the ToT calibration: ToT vs injected charge
	//(TOT = p0 (p1+Q)/(p2+Q))
	double Q = x[0];
	double tot = par[0] * (par[1] + Q) / (par[2] + Q);
	return tot;
}

ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::ToTCalib_FEI3_TurboDaq_converted_from_ToT_file() {
	calib = 0;
	correction = 1.0;
	calibration_type = "fei3_turbodaq_converted_from_tot_file";
}

ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::~ToTCalib_FEI3_TurboDaq_converted_from_ToT_file() {
	delete calib;
	if (hasPerPixelCalib()) {
		for (map<int, TF1*>::iterator i = calibpp.begin(); i != calibpp.end();
				++i) {
			if (i->second != NULL)
				delete i->second;
		}
	}
}

double ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::chargepp(const int& tot, const int& col, const int&row) {
	int i = index(col, row);
	if (calibpp.find(i) == calibpp.end()) {
		return charge(tot);
	} else {
		if (tot < calibpp[i]->GetMinimum()) {
			return calibpp[i]->GetX(calibpp[i]->GetMinimum());
		} else if (calibpp[i]->GetMaximum() < tot) {
			return calibpp[i]->GetX(calibpp[i]->GetMaximum());
		} else {
			return calibpp[i]->GetX(tot);
		}
	}
}

void ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::checkRange(const int& tot) {
	assert(calib != 0);
	double max = calib->GetMaximum();
	double min = calib->GetMinimum();
	if (min > tot || max < tot) {
		cout << "WARNING in totcalib: ToT " << tot << " is out of range ["
				<< calib->GetMinimum() << "," << calib->GetMaximum() << "]"
				<< endl;
		//exit(-1);
	}
}

double ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::charge(const int& tot) {
	//check the range
	assert(calib != NULL);
	checkRange(tot);
	if (tot < calib->GetMinimum()) {
		return calib->GetX(calib->GetMinimum());
	} else if (calib->GetMaximum() < tot) {
		return calib->GetX(calib->GetMaximum());
	} else {
		return calib->GetX(tot);
	}
}

bool ToTCalib_FEI3_TurboDaq_converted_from_ToT_file::addToTCalib(const std::string& filename, const int& iden,
		bool includePerPixel) {
	//  std::cout << "Adding ToT calibration from file " << filename << " to iden " << iden << endl;
	//Check that it doesn't exist
	if (not (calib == NULL)) {
		cerr
				<< "Error in totcalib: The calibration already exist. Function name: "
				<< calib->GetName() << endl;
		exit(-1);
	}

	TDirectory* dir = gDirectory;
	TFile* f = new TFile(filename.c_str());
	if (!f->IsOpen()) {
		cerr << "Error in totcalib: The file: " << filename
				<< " cannot be opened!" << std::endl;
		exit(-1);
	}

	//There are two calibration stored in some files
	//1. Calibrations averaging over all pixels
	//2. Calibration per pixel
	// Normally we use the averaged one.

	//Get the overall function
	TString name = TString::Format("ftotToQ_%i", iden);
	TF1* fnc = (TF1*) f->Get(name);
	if (!fnc) {
		cerr << "Error in totcalib: Function " << name
				<< " cannot be found in file:" << filename << std::endl;
		exit(-1);
	}
	//since it was fitted on range I have to make a new one to cover all possible values
	//this can probably be handled in some smarter way (tried setRange...) FIX THIS!
	calib = new TF1(TString::Format("ftotcalibave_%i", iden), conv, 0., 200000.,
			3);
	calib->SetParameters(fnc->GetParameter(0), fnc->GetParameter(1),
			fnc->GetParameter(2));
	cout << "Added ToT calibration function (" << calib->GetName() << ")"
			<< endl;
	if (includePerPixel) {
		cout << "Adding ToT calibration functions per pixel " << endl;
		//Get a the calibration function for each pixel
		int counter = 0;
		TF1* fnc2;
		for (int icol = 0; icol != cols; ++icol) {
			for (int irow = 0; irow != rows; ++irow) {
				TString s = TString::Format("ftotToQ_%i_%i_%i", irow, icol,
						iden);
				fnc = (TF1*) f->Get(s);
				int i = index(icol, irow);
				assert(calibpp.find(i) == calibpp.end());
				if (!fnc) {
					cout << "Warning in totcalib: pixel TOT calib for pixel "
							<< "row=" << irow << " col=" << icol
							<< " didn't exist!" << endl;
					assert(calibpp.find(i) == calibpp.end());
					calibpp[i] = NULL;
				} else {
					//same story with the range...
					fnc2 = new TF1(TString::Format("ftotcalibpp_%i", iden),
							conv, 0., 200000., 3);
					fnc2->SetParameters(fnc->GetParameter(0),
							fnc->GetParameter(1), fnc->GetParameter(2));
					calibpp[i] = fnc2;
					counter++;
				}
			}
		}
		cout << "ToTCalib: Added " << counter << " pixel calibration functions"
				<< std::endl;
	}
	f->Close();
	delete f;
	dir->cd();
	return true;
}
