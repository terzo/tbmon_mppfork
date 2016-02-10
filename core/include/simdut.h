#ifndef SIMDUT_H
#define SIMDUT_H

#include <string>
#include <vector>

#include "simPixelEdep.h"
#include "tbconfig.h"

/*! \brief SIMULATION: Interface/baseclass defining pixel digitization models.
 *
 */
class Simdut {
private:

public:
  //Stuff to be done before anything else
  virtual void init(TbConfig &config) {};
  //This function, which is implemented in the inheriting class,
  // makes digits for this event from the edeps.
  // The digits are stored in the vector hits, which is passed by reference.
  // This vector should be emptied *before* sending it to any model!
  virtual void digitize(Event& event,
                        const TbConfig& config,
                        std::vector<PllHit*>& hits) = 0;
  //Stuff to be done after all events have been processed
  virtual void finalize(const TbConfig &config) {};

  //Name of the model
  std::string modelName;
};


#endif //SIMDUT_H
