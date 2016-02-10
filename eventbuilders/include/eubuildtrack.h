#ifndef EUBUILDTRACK_H
#define EUBUILDTRACK_H

#include <map>
#include <algorithm>
#include <assert.h>

#include "TTree.h"

#include "event.h"
#include "Track.h"
#include "eventbuilder.h"

/*! \brief EuBuildTrack reads in data from tbtrack class (at the moment only zspix and eutracks branches)
 *
 *
 */
class EuBuildTrack: public EventBuilder {

private:
	/*! \brief Stores track shift per DUT and run generated by checkalign() and read in by translator()
	 *
	 */
	struct TransCalib {
		int firstRun;
		int lastRun;
		int iden;
		double shiftX;
		double shiftY;
	};

	std::vector<std::list<int> > matchWindows;
	map<int, std::vector<PllHit*> > m_hitMap;

	map<int, TransCalib> m_translations;
	list<TransCalib> m_calibs;
	list<int> m_idens;
	int m_nMatch;

	int pllCounter;

	vector<PllHit*> allPllHits;

	TTree* tracktree;

	//Track tree
	int t_nTrackParams;
	int t_euEv;
	std::vector<double> *t_posX;
	std::vector<double> *t_posY;
	std::vector<double> *t_dxdz;
	std::vector<double> *t_dydz;
	std::vector<int> *t_iden;
	std::vector<int> *t_trackNum;
	std::vector<double> *t_chi2;
	std::vector<double> *t_ndof;

	//Pixel tree
	TTree* pixeltree;

	int p_nHits;
	int p_euEv;
	std::vector<int> *p_col;
	std::vector<int> *p_row;
	std::vector<int> *p_tot;
	std::vector<int> *p_iden;
	std::vector<int> *p_lv1;
	std::vector<int> *p_chip;

	//Global rotation tree
	TTree* globalrotationtree;
	std::vector<int> *gr_ID;
	std::vector<double> *gr_Alpha;	/// Angle of sensor along x from alignment
	std::vector<double> *gr_Beta; 	/// Angle of sensor along y from alignment
	std::vector<double> *gr_Gamma;	/// Angle of sensor along z from alignment
	std::vector<double> *gr_RotXY;/// Angle of sensor in XY plane from GEAR file
	std::vector<double> *gr_RotZX;/// Angle of sensor in ZY plane from GEAR file
	std::vector<double> *gr_RotZY;/// Angle of sensor in ZX plane from GEAR file
	std::vector<double> *gr_RotXYErr;	/// Error on angle in XY
	std::vector<double> *gr_RotZXErr;	/// Error on angle in ZX
	std::vector<double> *gr_RotZYErr;	/// Error on angle in ZY

	bool hasRotationTree;

	//Store PLL's inorder to be able to delete them next go around
	std::vector<PllHit*> hits;

public:
	EuBuildTrack() {
		name = "EuBuildTrack";
	}

	int findMatches(Event &event);
	virtual void init(TbConfig &config);
	virtual void initRun(TbConfig &config);
	virtual void buildEvent(Event &event, map<int, Event> &events, TbConfig &);
	virtual void initEvent(TbConfig &config, map<int, Event> &events);
	void addTranslation(int iden, int firstRun, int lastRun, double shiftX,
			double shiftY);
	void addTranslation(const char* filename, int iden, double addX = 0,
			double addY = 0);
	void addMatchDUT(int iden);
	void nMatches(int nmatch);
	bool getHasRotationTree() {
		return hasRotationTree;
	}
	;
};

#endif //EUBUILDTRACK_H
