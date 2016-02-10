#ifndef SIMTRUTHHIT_H
#define SIMTRUTHHIT_H

#include <iostream>
#include <string>

#include "simThreeVector.h"

/*! \brief SIMULATION: Data storage for one simulation truth hit.
 *
 */
class simTruthHit {
public:
	simThreeVector pos;        //Global coordinate of particle
	simThreeVector posLocal; //Local -*- (transformed into DUT coordinate system)
	simThreeVector posLocalRaw;   //Local -*- (no transformations)
	simThreeVector momDir;     //Direction of particle momentum (global frame)
	simThreeVector momDirLocal;     // -*- (local frame)
	double kinE;               //Particle kinetic energy
	string planeID;            //Device name in TestBeamSim
	bool firstStep;          //Not used
	int particleID;         //Track number of particle in TestBeamSim
	int particleType;       //PDG ID of particle
	int stepNum;            //Step number for this track in TestBeamSim

	void print() {
		std::cout << "simTruthHit:" << std::endl;
		std::cout << "\t pos          =" << pos.toString() << std::endl;
		std::cout << "\t posLocal     =" << posLocal.toString() << std::endl;
		std::cout << "\t momDir       =" << momDir.toString() << std::endl;
		std::cout << "\t momDirLocal  =" << momDirLocal.toString() << std::endl;
		std::cout << "\t kinE         =" << kinE << std::endl;
		std::cout << "\t planeID      =" << planeID << std::endl;
		std::cout << "\t firstStep    =" << firstStep << std::endl;
		std::cout << "\t particleID   =" << particleID << std::endl;
		std::cout << "\t particleType =" << particleType << std::endl;
		std::cout << "\t stepNum      =" << stepNum << std::endl;
	}
};

#endif //SIMTRUTHHIT_H
