#ifndef CLUSTERS_H
#define CLUSTERS_H

#include "dut.h"
#include "event.h"
#include "tbutils.h"
#include "Track.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>

namespace cluster{
  const double pi = atan2(0.0f,-1.0f);

  // All clusters
  int    getSumTot( const event::cluster_t &cluster );
  int    getSumCharge( const event::cluster_t &cluster, const Event& event );
  int    getSumChargePP( const event::cluster_t &cluster, const Event& event);	
  //The following methods returns the index of a cluster
  int    getMaxTotCluster(const vector< event::cluster_t > &clusters);
  int    getMatched(const vector< event::cluster_t > &clusters, const Event &event);
  // Single cluster, returns index of cell
  int    getMaxTotCell(const event::cluster_t &clusters);

  enum ClusterPositionAlgorithm{
    kUnweighted,
    kChargeWeighted,
    kEtaCorrected,
    kMaxCell,
    kDigitalHeadTail1D,
    kAnalogHeadTail1D,
    kOneSidedAnalogHeadTail1D,
    kUnknown //Make sure that kUnknown stays last
  };

  static const int NumClusterPositionAlgorithm = 8;
  
  /// Names of cluster position finding algorithms. 
  /// Words in square brackets will be used as histo names. The rest as title
  static const string ClusterPositionAlgorithmNames[NumClusterPositionAlgorithm]={
    "[geomean]Unweighted", //Digital
    "[qmean]Charge Weighted", //Analog
    "[eta]#eta Corrected",
    "[maxtot]Max ToT",
    "[DHT]Digital Head Tail1 1D",
    "[AHT]Analog Head Tail 1D",
    "[1SAHT]One Sided Analog Head Tail 1D",
    "[Unknown]Unknown"};

    /*
  static std::map<std::string,int> ClusterPositionAlgorithmShortNameMapping = {
		  {"kUnweighted",kUnweighted},
		  {"kChargeWeighted",kChargeWeighted},
		  {"kEtaCorrected",kEtaCorrected},
		  {"kMaxCell",kMaxCell},
		  {"kDigitalHeadTail1D",kDigitalHeadTail1D},
		  {"kAnalogHeadTail1D",kAnalogHeadTail1D},
		  {"kOneSidedAnalogHeadTail1D",kOneSidedAnalogHeadTail1D},
		  {"kUnknown",kUnknown}
  };
    */

  static int lineColor[cluster::NumClusterPositionAlgorithm] = {kBlue, kRed, 
    kBlack, kOrange, kGreen, kMagenta, kViolet, kSpring};
  static int areaColor[cluster::NumClusterPositionAlgorithm] = {kWhite, kWhite, 
    kBlack, kOrange, kWhite, kWhite, kWhite, kWhite};
  
  void getAlgNameTitle(int alg, string* name, string *title = 0);
  string getAlgName(int alg);
  string getAlgTitle(int alg);
  int shortAlgoNameToID(std::string input);

  /*! \brief get optimal cluster center finding algorithm
   *
   * get optimal cluster center finding algorithm according to:
   * R. Turchetta, Spatial resolution of silicon microstrip detectors,
   * Nuclear Instruments and Methods in Physics Research A 335 (1993) 44-58
   */
  int getOptimalClusterCenterAlgorithm(double angle);
   
  double getCol(const event::cluster_t &cluster, const Event& event, int algo);
  double getRow(const event::cluster_t &cluster, const Event& event, int algo);
  double getX(const event::cluster_t &cluster, const Event& event, int algo);
  double getY(const event::cluster_t &cluster, const Event& event, int algo);
 

  //********One dimensional position finding algorithms**************//
  //Cols and rows
  double getUnWeightedCol(const event::cluster_t &cluster);
  double getUnWeightedRow(const event::cluster_t &cluster);
  //center-of-gravity
  double getChargeWeightedCol(const event::cluster_t &cluster);
  double getChargeWeightedRow(const event::cluster_t &cluster);
  double getEtaCorrectedCol(const event::cluster_t &cluster, DUT* dut);
  double getEtaCorrectedRow(const event::cluster_t &cluster, DUT* dut);
  // cell with maximum ToT
  double getMaxTotCellCol(const event::cluster_t &cluster);
  double getMaxTotCellRow(const event::cluster_t &cluster);
  //digital head-tail (DHT) PFA
  double getClusterCenterColByDigitalHeadTail1D(const event::cluster_t &cluster);
  double getClusterCenterRowByDigitalHeadTail1D(const event::cluster_t &cluster);
  //analog head-tail (AHT) PFA
  double getClusterCenterColByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);
  double getClusterCenterRowByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);
  //one-sided analog head-tail (1S-AHT) PFA
  double getClusterCenterColByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);
  double getClusterCenterRowByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);


  // X and Y in µm
  double getUnWeightedX( const event::cluster_t &cluster, const Event& event );
  double getUnWeightedY( const event::cluster_t &cluster, const Event& event );
  //center-of-gravity
  double getChargeWeightedX( const event::cluster_t &cluster, const Event& event );
  double getChargeWeightedY( const event::cluster_t &cluster, const Event& event );
  double getEtaCorrectedX(const event::cluster_t &cluster, const Event& event, DUT* dut);
  double getEtaCorrectedY(const event::cluster_t &cluster, const Event& event, DUT* dut);
  // cell with maximum ToT
  double getMaxTotCellX(const event::cluster_t &cluster, const Event& event);
  double getMaxTotCellY(const event::cluster_t &cluster, const Event& event);  
  //digital head-tail (DHT) PFA
  double getClusterCenterXByDigitalHeadTail1D(const event::cluster_t &cluster, const Event& event );
  double getClusterCenterYByDigitalHeadTail1D(const event::cluster_t &cluster, const Event& event );
  //analog head-tail (AHT) PFA
  double getClusterCenterXByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);
  double getClusterCenterYByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);
  //one-sided analog head-tail (1S-AHT) PFA
  double getClusterCenterXByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);
  double getClusterCenterYByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event);

  //*****************************************************************//

  //Does track match pixel hit?
  bool isMatch(const event::cluster_t &cluster, const Event &event, int algo=kUnknown);
  bool isCentralRegion(const event::cluster_t &cluster, const Event &event);
  int spannedRows(const event::cluster_t &cluster);
  int spannedCols(const event::cluster_t &cluster);
  
};

#endif //CLUSTERS_H
