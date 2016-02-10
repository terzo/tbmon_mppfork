#ifndef CHECKREGION_H
#define CHECKREGION_H

// standard header files
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief Checks whether track is in the central region, requires addMask()
 *
 */
class CheckRegion : public EventBuilder {
private:
	double maskRad;	///< Radius in units of pixel pitch that will be employed for masking tracks in the vicinity of masked pixels

public:
  CheckRegion () {
    name = "CheckRegion";
  }
  virtual void init(TbConfig &config);
  virtual void initRun(TbConfig &config){;};
  virtual void buildEvent(Event &event, map<int,Event> &events, TbConfig &);
  bool BadRegion(Event &event, TbConfig &config);
};

#endif //CHECKREGION_H
