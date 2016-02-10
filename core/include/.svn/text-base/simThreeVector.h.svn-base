#ifndef SIMTHREEVECTOR_H
#define SIMTHREEVECTOR_H

#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

class Event;

using namespace std;

//Forward declarations needed b.c. strange syntax
class simThreeVector;
istream& operator>>(istream& is, simThreeVector& vec) throw (ios_base::failure);

/*! \brief SIMULATION: Used to parse three-vectors from text files, format (x,y,z)
 *
 * Kyrre N. Sjobak
 */
class simThreeVector {
public:
  double data [3];

  friend istream& operator>>(istream& is, simThreeVector& vec) throw (ios_base::failure);

  string toString();

  void localTranslate(const Event& event);
};


#endif //SIMTHREEVECTOR_H
