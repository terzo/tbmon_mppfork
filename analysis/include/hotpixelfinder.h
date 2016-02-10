#ifndef HOTPIXELFINDER_H
#define HOTPIXELFINDER_H

//tbmon header files
#include "event.h"
#include "tbanalysis.h"

//root header files
#include "TH2D.h"
#include "TExec.h"

//standard header files
#include <fstream>
#include <iostream>
#include <stdlib.h>

/*! \brief HotPixelFinder is an analysis class which creates mask files for hot and dead pixels
 *
 * HotPixelFinder is an analysis class which creates mask files thus masking pixels
 * which do not show activity at all (called dead) or which are firing too often (above a certain
 * occupancy threshold). The occupancy threshold can be changed using the command line parameter
 * -P:A_hotpixelfinder_maximumOccupancy 0.1
 * Occupancy is the ratio of hits outside the lv1 window which is defined in the driver.cc and all hits.
 * The analysis creates one global mask file per DUT.
 */
class HotPixelFinder: public TbAnalysis {

private:
	TH2D* h_hitmap; ///< all hits of one DUT (of events where the track was accepted by the referencing) are filled into this 2D histogramm to create a hit map
	TH2D* h_masked; ///< all masked pixels, set in the two helper functions checkDead() and checkNoise(), 0=active, 1=masked as dead, 2=masked as noisy
	TH2D* h_outOfTime; ///< map with all hits which are out of the designated time window
	TH1D* h_lv1; ///< level 1 timing distribution of one DUT
	TH1D* h_noiseOccupancy; ///< noise occupancy distribution (filled from the out of time histogram wherein it is calculated in checkNoise()

	int eventCount; ///< total number of events
	int lv1Min, lv1Max; ///< level 1 timing window, within hits are "in time"
	int nCols; ///< number of columns in sensor - gets filled from the TBConfig object
	int nRows; ///< number of rows in sensor - gets filled from the TBConfig object

	double _maximumOccupancy; ///< minimum occupancy for flagging a pixel as "hot"

	int numberOfHotPixels; ///< total number of hot pixels in one DUT
	int numberOfDeadPixels; ///< total number of dead pixels in one DUT

	ofstream stream; ///< file handle used to (re-)create and write the mask file in text format

	void checkDead(); ///< helper function which looks at the filled hitmap histogram and searches for non-active pixels
	void checkNoise(); ///< helper function which looks at the filled out of time histogram, calculates an occupancy for each pixels and applies the cut

public:
	virtual void init(TbConfig &config);
	virtual void event(const TbConfig &config, const Event &event);
	virtual void finalize(const TbConfig &config);
};

#endif //HOTPIXELFINDER_H
