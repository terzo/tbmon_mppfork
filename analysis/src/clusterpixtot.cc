#include <clusterpixtot.h>
#include <clusters.h>

#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

using namespace std;

void ClusterPixToT::init(TbConfig &config) {
	// Add command line option to restrict cluster size for histos
	string key;
	key = Form("clupixtot_%u_Xsize", this->iden);
	_fixedClusterSizeX = config.cmdLineExtras_argGetter(key, (double) 5,
			"Restrict cluster X size for histograms(default: 5)");

//    key = Form("clupixtot_%u_Ysize", this->iden);
//    _fixedClusterSizeY = config.cmdLineExtras_argGetter(key, (double)1,
//          "Restrict cluster Y size for histograms(default: 1)");

}

void ClusterPixToT::initRun(const TbConfig &config) {
	string streamName(config.getOutStreamName(name, "clustertottree"));
	streamName = streamName.substr(0, streamName.find("."));
	streamName = streamName + ".root";

	_outTree = new TTree(Form("ToTtree_%d", iden), Form("ToTtree_%d", iden));
	_outTree->Branch("clusterNumber", &_clusterNumber);
	_outTree->Branch("pixelToT", &_pixelToT);
	_outTree->Branch("pixelPosX", &_pixelPosX);
	_outTree->Branch("pixelPosY", &_pixelPosY);
	_outTree->Branch("clusterSizeX", &_clusterSizeX, "clusterSizeX/I");
	_outTree->Branch("clusterSizeY", &_clusterSizeY);

	_1stPixToTdependentToT.clear();
	_clusterNumber = -1;
}

void ClusterPixToT::event(const TbConfig &config, const Event &event) {
	// Check that track has passed all cuts
	if (event.fTrack != event::kGood) {
		return;
	}
	// Check that clusters have been made successfully
	if (event.fClusters != event::kGood) {
		return;
	}
	// Are we in a good region of the chip?
	if (event.fTrackRegion != event::kGood) {
		return;
	}
	// Get cluster matched to track, if not found index is -1
	int matched = cluster::getMatched(event.clusters, event);
	if (matched == -1) {
		return;
	}

	int minPosX = event.dut->getNcols();
	int maxPosX = 0;
	int minPosY = event.dut->getNrows();
	int maxPosY = 0;
	int minPosToT = -1;

	// Calculate cluster boundary and save ToT of pixel with minimal X position
	for (vector<PllHit*>::const_iterator jj =
			event.clusters.at(matched).begin();
			jj != event.clusters.at(matched).end(); ++jj) {
		if ((*jj)->col > maxPosX) {
			maxPosX = (*jj)->col;
		}
		if ((*jj)->col < minPosX) {
			minPosX = (*jj)->col;
			minPosToT = (*jj)->tot;
		}

		if ((*jj)->row > maxPosY) {
			maxPosY = (*jj)->row;
		}
		if ((*jj)->row < minPosY) {
			minPosY = (*jj)->row;
		}
	}

	_clusterSizeY = maxPosY - minPosY + 1;
	_clusterSizeX = maxPosX - minPosX + 1;
	_clusterNumber++;

	// Save position and ToT of every pixel in cluster
	for (vector<PllHit*>::const_iterator jj =
			event.clusters.at(matched).begin();
			jj != event.clusters.at(matched).end(); ++jj) {
		_pixelToT = (*jj)->tot;
		_pixelPosX = (*jj)->col - minPosX;
		_pixelPosY = (*jj)->row - minPosY;
		_outTree->Fill();

		// If cluster size passes cut add pixel to ToT collection
		if (_clusterSizeY == _fixedClusterSizeY
				&& _clusterSizeX == _fixedClusterSizeX) {
			if (_pixelToT < 14) {
				_1stPixToTdependentToT[minPosToT][_pixelPosX][_pixelToT]++;
			}

			// Book _singlePixToTDistr histo for pixel position if necessary and fill it
			TH1F *singlePixHist = _singlePixToTDistr[_pixelPosX];
			if (!singlePixHist) {
				TBALOG(kINFO) << "Booking histogram "
						<< Form("totdistr_dut%d_pix_%d", iden, _pixelPosX)
						<< endl;
				_singlePixToTDistr[_pixelPosX] = singlePixHist = new TH1F(
						Form("totdistr_dut%d_pix_%d", iden, _pixelPosX),
						Form("ToT distribution for pixel %d (DUT %d)",
								_pixelPosX, iden), 15, -0.5, 14.5);
				singlePixHist->GetXaxis()->SetTitle("Charge [ToT]");
			}
			singlePixHist->Fill(_pixelToT);
		}
	}
}

void ClusterPixToT::finalizeRun(const TbConfig &config) {
	// Iterate over every ToT of the first pixel (minimum X position) in cluster
	for (map<int, map<int, map<int, int> > >::iterator pix1ToTiter =
			_1stPixToTdependentToT.begin();
			pix1ToTiter != _1stPixToTdependentToT.end(); pix1ToTiter++) {
		int minPosToT = pix1ToTiter->first;

		//Iterate over every position of pixels in cluster and calculate their mean ToT
		for (map<int, map<int, int> >::iterator pixPosIter =
				pix1ToTiter->second.begin();
				pixPosIter != pix1ToTiter->second.end(); pixPosIter++) {
			int pixNum = pixPosIter->first;

			double mean = 0;
			double sigma = 0;
			int N = 0;
			for (map<int, int>::iterator ToTiter = pixPosIter->second.begin();
					ToTiter != pixPosIter->second.end(); ToTiter++) {
				mean += ToTiter->first * ToTiter->second;
				N += ToTiter->second;
			}
			if (N != 0) {
				mean /= N;
			} else {
				TBALOG(kERROR) << "Number of elements == 0" << endl;
				return;
			}

			for (map<int, int>::iterator totit = pixPosIter->second.begin();
					totit != pixPosIter->second.end(); totit++) {
				sigma += totit->second * pow(totit->first - mean, 2);
			}
			if (N > 1) {
				sigma = sqrt(1. / (N - 1) * sigma);
			} else {
				sigma = 0;
			}

			TBALOG(kDEBUG) << "1st pix tot: " << minPosToT << " pixPos: "
					<< pixNum << ": " << mean << "+-" << sigma << " (" << N
					<< ")" << endl;

			if (!_1stPixDepTotDistr[minPosToT]) { //C++ Standard, 8.5 paragraph 5
				TBALOG(kINFO) << "Booking histogram "
						<< Form("totdistr_dut%d_pix1tot_%d", iden, minPosToT)
						<< endl;
				_1stPixDepTotDistr[minPosToT] = new TH1F(
						Form("totdistr_dut%d_pix1tot_%d", iden, minPosToT),
						Form("ToT distribution DUT %d (1st pixel ToT=%d)", iden,
								minPosToT), 20, -.5, 19.5);
				_1stPixDepTotDistr[minPosToT]->GetXaxis()->SetTitle(
						"Pixel position [pitch]");
				_1stPixDepTotDistr[minPosToT]->GetYaxis()->SetTitle(
						"Charge [ToT]");
				_1stPixDepTotDistr[minPosToT]->GetYaxis()->SetTitleOffset(0.8);
			}
			int bin = _1stPixDepTotDistr[minPosToT]->Fill(pixNum, mean);
			_1stPixDepTotDistr[minPosToT]->SetBinError(bin, sigma);

			if (!_1stPixDepSinglePixTotDistr[pixNum]) {
				TBALOG(kINFO) << "Booking histogram "
						<< Form("totdistr_dut%d_pixpos_%d", iden, pixNum)
						<< endl;
				_1stPixDepSinglePixTotDistr[pixNum] = new TH1F(
						Form("1st_pix_dep_totdistr_dut%d_pix_%d", iden, pixNum),
						Form("1st pixel dependent ToT (pix %d, DUT %d)", pixNum,
								iden), 15, -0.5, 14.5);
				_1stPixDepSinglePixTotDistr[pixNum]->GetXaxis()->SetTitle(
						"Charge of first pixel [ToT]");
				_1stPixDepSinglePixTotDistr[pixNum]->GetYaxis()->SetTitle(
						Form("Charge of pixel %d[ToT]", pixNum));
				_1stPixDepSinglePixTotDistr[pixNum]->GetYaxis()->SetTitleOffset(
						0.8);
			}
			bin = _1stPixDepSinglePixTotDistr[pixNum]->Fill(minPosToT, mean);
			_1stPixDepSinglePixTotDistr[pixNum]->SetBinError(bin, sigma);
		}
	}

	// Save histos and tree
	config.saveToFile(this->name, _outTree->GetName(), _outTree);

	for (map<int, TH1F*>::iterator it = _singlePixToTDistr.begin();
			it != _singlePixToTDistr.end(); it++) {
		config.saveToFile(this->name, it->second->GetName(), it->second);
	}

	for (map<int, TH1F*>::iterator it = _1stPixDepTotDistr.begin();
			it != _1stPixDepTotDistr.end(); it++) {
		config.saveToFile(this->name, it->second->GetName(), it->second);
	}

	for (map<int, TH1F*>::iterator it = _1stPixDepSinglePixTotDistr.begin();
			it != _1stPixDepSinglePixTotDistr.end(); it++) {
		config.saveToFile(this->name, it->second->GetName(), it->second);
	}
}
