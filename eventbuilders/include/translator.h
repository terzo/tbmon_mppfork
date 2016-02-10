#ifndef TRANSLATOR_H
#define TRANSLATOR_H

// standard header files
#include <map>

// tbmon header files
#include "event.h"
#include "eventbuilder.h"
#include "Track.h"

/*! \brief Reads shift per DUT and run from file (generated by checkalign) and applies it to each track.
 *
 */
class Translator: public EventBuilder {
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
	map<int, TransCalib> m_translations;
	list<TransCalib> m_calibs;
public:
	Translator() {
		name = "Translator";
	}
	;
	/** @fn load in translations stored for the individual DUTs */
	virtual void initRun(TbConfig &config);
	/** @fn translate event.trackX/trackY by the translation belonging to the DUT of the event in question */
	virtual void buildEvent(Event &event, map<int, Event> &events, TbConfig &);
	void addTranslation(int iden, int firstRun, int lastRun, double shiftX,
			double shiftY);
	void addTranslation(const char* filename, int iden, double addX = 0,
			double addY = 0);
};

#endif //TRANSLATOR_H