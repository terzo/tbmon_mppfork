#ifndef CHECKALIGN_RUNNINGXY_H
#define CHECKALIGN_RUNNINGXY_H

// standard header files
#include <vector>
#include <string>
#include <iostream>
#include <deque>

// root header files
#include <TProfile.h>
#include <TFile.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TTree.h>
#include <TROOT.h>

// tbmon header files
#include "clusters.h"
#include "event.h"
#include "tbanalysis.h"

/*! \brief Calculates residuals for a sliding window of events and generates fromt hose residuals a
 * time dependent translation correction which can be read in using translationXY
 *
 */
class CheckAlignRunningXY: public TbAnalysis {

private:

	/*! \brief Helper class for CheckAlignRunningXY(), see TbTupleAna
	 *
	 */
	class SensorInfoObject {
	protected:
		int _sensorId;

	public:
		SensorInfoObject() {
		}
		;
		SensorInfoObject(const int sensorId) :
				_sensorId(sensorId) {
		}
		;
		virtual ~SensorInfoObject() {
		}
		;
		int getId() {
			return _sensorId;
		}
		virtual void clear() {
		}
	};

	/*! \brief Helper class for CheckAlignRunningXY(), see TbTupleAna
	 *
	 */
	class SensorInfoObjectCollection {
	protected:
		std::map<int, SensorInfoObject*> _infoMap;

	public:
		SensorInfoObjectCollection() {
		}
		;

		void add(SensorInfoObject *obj) {
			std::cout << "Adding SensorInfo with Id " << obj->getId() << std::endl;
			_infoMap[obj->getId()] = obj;
		}

		SensorInfoObject* get(const int id) {
			std::map<int, SensorInfoObject*>::iterator it = _infoMap.find(id);
			if (it == _infoMap.end())
				return NULL;
			return it->second;
		}

		std::vector<SensorInfoObject*> get() {
			std::map<int, SensorInfoObject*>::iterator it;
			std::vector<SensorInfoObject*> out;
			for (it = _infoMap.begin(); it != _infoMap.end(); it++) {
				out.push_back(it->second);
			}
			return out;
		}
		virtual void clear() {
			std::map<int, SensorInfoObject*>::iterator it = _infoMap.begin();
			for (it; it != _infoMap.end(); it++) {
				it->second->clear();
				delete it->second;
			}
			_infoMap.clear();
		}

		virtual ~SensorInfoObjectCollection() {
			clear();
		}

	};

	/*! \brief Helper class for CheckAlignRunningXY(), see TbTupleAna
	 *
	 */
	class SensorInfo: public SensorInfoObject {
	public:
		enum SensorType {
			kUnknown = 0, kM26 = 1, kAPIXFeI3 = 2, kAPIXFeI4 = 3
		};
	protected:

		SensorType _type;
		std::string _prefix;
		std::string _title;
		int _pixX;
		int _pixY;
		float _pitchX;
		float _pitchY;


		int getPixelSerialized(const int x, const int y) {
			return x + y * _pixX;
		}

		void getPixelSerialized(const int i, int &x, int &y) {
			y = i / _pixX;
			x = i % _pixX;
		}

	public:

		//! Copy  constructor
		SensorInfo(const SensorInfo &s) {
			_type = s._type;
			_prefix = s._prefix;
			_title = s._title;
			_sensorId = s._sensorId;
			_pixX = s._pixX;
			_pixY = s._pixY;
			_pitchX = s._pitchX;
			_pitchY = s._pitchY;
		}
		SensorInfo() :
				_type(kUnknown) {
		}
		;
		SensorInfo(const int id) :
				SensorInfoObject(id), _type(kUnknown) {
		}
		SensorInfo(const int id, const SensorType type,
				const std::string title = "", const std::string prefix = "");

		virtual ~SensorInfo() {
		}
		void getPitch(float &pitchX, float &pitchY) {
			pitchX = _pitchX;
			pitchY = _pitchY;
		}
		void getNPixel(int &pixX, int &pixY) {
			pixX = _pixX;
			pixY = _pixY;
		}
		std::string getPrefix() {
			return _prefix;
		}
		std::string getTitle() {
			return _title;
		}
		SensorType getType() {
			return _type;
		}
	};

	/*! \brief Helper class for CheckAlignRunningXY(), see TbTupleAna
	 *
	 */
	class RunningAverage: public SensorInfoObject {
	protected:
		//! The number of entries to calculate the mean from
		int _width;
		//! The queue containing the last entries
		std::deque<double> _lastValues;
		//! This is the sum of all values in the queue
		double _sum;
		//! The event at which the queue was filled first
		int _meanBeginEvent;

	public:
		//! Constructor
		RunningAverage(const int sensorId, const int width) :
				SensorInfoObject(sensorId), _width(width), _sum(0), _meanBeginEvent(
						-1) {
		}
		//! Add a value to the queue
		/**
		 * if the queue is filled then delete the first element, so that the number of elements in the queue is always constant
		 **/
		void add(double value) {
			if (_lastValues.size() < _width) {
				_sum = _sum + value;
				_lastValues.push_back(value);
			} else {
				_sum = _sum - _lastValues.front();
				_lastValues.pop_front();
				_sum = _sum + value;
				_lastValues.push_back(value);
			}
		}
		//! Checks if the meanBegin variable is set
		bool isMeanBeginSet() {
			return (_meanBeginEvent != -1);
		}
		//! Setter for the MeanBeginVariable
		void setMeanBegin(const int meanBegin) {
			_meanBeginEvent = meanBegin;
		}
		//! Getter for the meanBeginVariable
		/**
		 * @ return -1 if meanBegin is not set. Check this with isMeanBeginSet()
		 * @return meanBegin variable otherwise
		 **/
		int getMeanBegin() {
			return _meanBeginEvent;
		}
		//! Calculated Mean value
		/**
		 * this function will return the mean in anycase. Check if the queue is filled with the isFilledUp() function.
		 **/
		double getMean() {
			if (_lastValues.size() == 0)
				return 0;
			return _sum / _lastValues.size();

		}

		double getFillRatio() {
			return (double) _lastValues.size() / _width;
		}
		//! Checks if the queue is filled up
		bool isFilledUp() {
			return !(_lastValues.size() < _width);
		}
		//! Calculation of the first n Events
		/**
		 * This method only works is at this moment the queue is completely filled first time (When the queue is in the state that setMeanBegin() is just called)
		 **/
		double getMeanBegin(const int index) {
			if (!isFilledUp())
				return 0.;
			std::deque<double>::iterator it;
			int myIndex = 0;
			int end = _lastValues.size();
			float ratioEvTrk = (float) _meanBeginEvent / end;

			double ret = 0;
			it = _lastValues.begin();
			while (myIndex < end) {
				//for (it = _lastValues.begin(); it != _lastValues.end(); it++) {
				if (myIndex < (int) (index / ratioEvTrk)) {
					ret += 2 * (*it);
					myIndex += 2;
				} else {
					ret += (*it);
					myIndex++;
				}
				it++;
			}

			return ret / end;

		}

	};

	/*! \brief Helper class for CheckAlignRunningXY(), see TbTupleAna
	 *
	 */
	class PostShiftAlign {
	public:
		//! The tree to store to or read from
		TTree *_tree;
		int _sensorId;
		float _shiftX;
		float _shiftY;
		int _eventNum;
		int _subsetSize;

		TBranch *_bSensorId;
		TBranch *_bShiftX;
		TBranch *_bShiftY;
		TBranch *_bEventNum;
		TBranch *_bSubsetSize;

		//!Constructor to write out the database into a TTree
		PostShiftAlign() :
				_sensorId(-1), _shiftX(0), _shiftY(0), _eventNum(0), _subsetSize(0) {
			_tree = new TTree("PostShiftAlign", "PostShiftAlign");
			//if (file) _tree->SetDirectory(file);
			_tree->Branch("sensorId", &_sensorId);
			_tree->Branch("shiftX", &_shiftX);
			_tree->Branch("shiftY", &_shiftY);
			_tree->Branch("eventNum", &_eventNum);
			_tree->Branch("subsetSize", &_subsetSize);
		}
		//! Constructor to read in a database from a TTree in a TFile
		PostShiftAlign(TFile *file) {
			_tree = (TTree*) file->Get("PostShiftAlign");
			if (!_tree)
				std::cerr << "Tree PostShiftAlign not found!" << std::endl;
			_tree->SetBranchAddress("sensorId", &_sensorId, &_bSensorId);
			_tree->SetBranchAddress("shiftX", &_shiftX, &_bShiftX);
			_tree->SetBranchAddress("shiftY", &_shiftY, &_bShiftY);
			_tree->SetBranchAddress("eventNum", &_eventNum, &_bEventNum);
			_tree->SetBranchAddress("subsetSize", &_subsetSize, &_bSubsetSize);

		}
		//! Set the directory of the TTree, so that Write() writes at the correct place
		void setDirectory(TFile * file) {
			_tree->SetDirectory(file);
		}
		//! Get the number of entries
		int getN() {
			return _tree->GetEntriesFast();
		}

		//! Load on entry into the class' member variables
		Int_t getEntry(Long64_t entry) {
			return _tree->GetEntry(entry);
		}

		//! Add the arguments to the TTree - This is used to fill the TTree
		void add(const int sensor, const float x, const float y, const int eventNum,
				const int size) {
			_sensorId = sensor;
			_shiftX = x;
			_shiftY = y;
			_eventNum = eventNum;
			_subsetSize = size;
			_tree->Fill();
		}

		//! Write out the TTree into a TFile
		void write(TFile *file = NULL) {
			if (file)
				_tree->SetDirectory(file);
			_tree->Write();

		}

	};

	TFile *_postShiftFile;

	int _postShiftSubset;
	int _firstPostShiftEvent;
	int _postShiftEventSkip;

	PostShiftAlign* _postShift;

	SensorInfoObjectCollection* _runningMeanX;
	SensorInfoObjectCollection* _runningMeanY;

	RunningAverage *avX;
	RunningAverage *avY;

	bool doCuts;

public:
	virtual void init(TbConfig &config);
	virtual void initRun(const TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalizeRun(const TbConfig &config);
	virtual void finalize(const TbConfig &config);
};

#endif //CHECKALIGN_RUNNINGXY_H
