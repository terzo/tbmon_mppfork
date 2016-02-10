#ifndef DUTMAKER_H
#define DUTMAKER_H

#include "dut.h"

DUT* makeFEI3(const char* name, int iden, int numElec, double epitchX=600.0, double epitchY=50.0, double refLimitX = 1.5, double refLimitY = 1.5, double thickness=250e-6){
  return( new DUT(name, iden, numElec, 400.0, 50.0, epitchX, epitchY, 18, 160, 17, 159, 3, 16, /*400 */ 800.0 , 150.0 /*50.0*/ /*100*/, refLimitX, refLimitY, thickness));
}

DUT* makeFEI4(const char* name, int iden, int numElec, double epitchX=250.0, double epitchY=50.0, double refLimitX = 1.5, double refLimitY = 1.5, double thickness=250e-6){ 
  return( new DUT(name, iden, numElec, 250.0, 50.0, epitchX, epitchY, 80, 336, 79, 335, 3, 16, /*500*/400.0, 150.0, refLimitX, refLimitY, thickness));
  // 400.0, 150.0
}

DUT* makeFEI4_25(const char* name, int iden, int numElec, double epitchX=500.0, double epitchY=25.0, double refLimitX = 1.5, double refLimitY = 1.5, double thickness=280e-6){ 
  return( new DUT(name, iden, numElec, 500.0, 25.0, epitchX, epitchY, 40, 672, 39, 671, 3, 16, 1000.0, 75.0, refLimitX, refLimitY, thickness));
}

#endif //DUTMAKER_H
