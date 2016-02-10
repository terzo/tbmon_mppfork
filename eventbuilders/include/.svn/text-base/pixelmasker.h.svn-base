#ifndef PIXELMASKER_H
#define PIXELMASKER_H

// standard header files
#include <cmath>
#include <vector>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief Applies pixel mask saved in DUT object that was previously read in by AddMasks() method.
 *
 * Applies pixel mask saved in DUT object that was previously read in by AddMasks() method. Not masked hits
 * are copied from rawhits to hits if they are within the previously defined LV1 timing window.
 */
class PixelMasker: public EventBuilder {
public:
	PixelMasker() {
		name = "PixelMasker";
	}
	virtual void buildEvent(Event &event, map<int, Event>&, TbConfig &);
};

#endif //PIXELMASKER_H
