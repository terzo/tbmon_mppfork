#ifndef SIMPIXELEDEP_H
#define SIMPIXELEDEP_H

#include <iostream>

#include "simThreeVector.h"

/*! \brief SIMULATION: Data storage for one simulation edep.
 *
 */
class simPixelEdep {
public:
  double edep;
  simThreeVector pos;
  simThreeVector posLocal;

  void print() {
    std::cout << "simPixelEdep:"                      << std::endl;
    std::cout << "\t edep   =" << edep                << std::endl;
    std::cout << "\t pos    =" << pos.toString()      << std::endl;
    std::cout << "\t edep   =" << posLocal.toString() << std::endl;
  }
};

#endif //SIMPIXELEDEP_H
