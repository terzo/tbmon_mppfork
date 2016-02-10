#include "clusters.h"

void cluster::getAlgNameTitle(int alg, string* name, string *title){
  if (alg < 0 || alg > kUnknown) {
    return;
  }
  
  bool nameFound = false;
  int endIndex = -1;
  string algName = ClusterPositionAlgorithmNames[alg];
  if (algName[0] == '['){
    endIndex = algName.find("]");
    if (endIndex != string::npos){
      nameFound = true;
    }
  }
  
  if (nameFound){
    if (name){
      *name = algName.substr(1, endIndex-1);
    }
    if (title){
      *title = algName.substr(endIndex+1, string::npos);
    }
  }else{
    if (title){
      *title = algName;
    }
    if (name){
      algName.erase(remove_if(algName.begin(), algName.end(), ::isspace), algName.end());
      *name = algName;
    }
  }
};

string cluster::getAlgName(int alg){
  string ret;
  getAlgNameTitle(alg, &ret, 0);
  return ret;
}

string cluster::getAlgTitle(int alg){
  string ret;
  getAlgNameTitle(alg, 0, &ret);
  return ret;
}

int cluster::getOptimalClusterCenterAlgorithm(double angle){
	//FIXME: stub
	std::cout << "error: this is a stub" << std::endl;
}

int cluster::shortAlgoNameToID(std::string input){
	std::map<std::string,int> ClusterPositionAlgorithmShortNameMapping;
	ClusterPositionAlgorithmShortNameMapping["kUnweighted"] = kUnweighted;
	ClusterPositionAlgorithmShortNameMapping["kChargeWeighted"] = kChargeWeighted;
	ClusterPositionAlgorithmShortNameMapping["kEtaCorrected"] = kEtaCorrected;
	ClusterPositionAlgorithmShortNameMapping["kMaxCell"] = kMaxCell;
	ClusterPositionAlgorithmShortNameMapping["kDigitalHeadTail1D"] = kDigitalHeadTail1D;
	ClusterPositionAlgorithmShortNameMapping["kAnalogHeadTail1D"] = kAnalogHeadTail1D;
	ClusterPositionAlgorithmShortNameMapping["kOneSidedAnalogHeadTail1D"] = kOneSidedAnalogHeadTail1D;
	ClusterPositionAlgorithmShortNameMapping["kUnknown"] = kUnknown;
		  
	if(ClusterPositionAlgorithmShortNameMapping.count(input)>0){
		int result =  ClusterPositionAlgorithmShortNameMapping[input];
		return result;
	} else
		return ClusterPositionAlgorithmShortNameMapping["kUnknown"];
}

int cluster::getSumTot( const event::cluster_t &cluster){
  int sumTot(0);
  for(int ii = 0; ii < cluster.size(); ii++){
    sumTot += cluster.at(ii)->tot;
  }
  return(sumTot);
}
int cluster::getSumCharge( const event::cluster_t &cluster, const Event& event){
  double sumCharge = 0.0;
  for(int ii = 0; ii < cluster.size(); ii++){
    sumCharge += event.dut->q(cluster.at(ii)->tot);
  }
  return((int) sumCharge);
}
int cluster::getSumChargePP( const event::cluster_t &cluster, const Event& event){
  double sumCharge = 0.0;
  for(int ii = 0; ii < cluster.size(); ii++){
    sumCharge += event.dut->q(cluster.at(ii)->tot,cluster.at(ii)->col,cluster.at(ii)->row);
  }
  return((int) sumCharge);
}
int cluster::getMaxTotCluster(const vector< event::cluster_t > &clusters){
  int index(-1), maxTot(0);
  for(int ii = 0; ii < clusters.size(); ii++){
    int tmpTot = getSumTot( clusters.at(ii) );
    if(tmpTot < maxTot) { continue; }
    index = ii;
    maxTot = tmpTot;
  }
  return(index);
}

int cluster::getMaxTotCell(const event::cluster_t &cluster){
  int index(-1), maxTot(0);
  for(int ii = 0; ii < cluster.size(); ii++){
    int tmpTot = cluster.at(ii)->tot;
    if(tmpTot < maxTot){ continue;}
    index = ii;
    maxTot = tmpTot;
  }
  return(index);
}

double cluster::getUnWeightedCol(const event::cluster_t &cluster){
  int n(0); double pos(0);
  for(int ii = 0; ii < cluster.size(); ii++){
    pos += cluster.at(ii)->col; n++;
  }
  return(pos/double(n));
}
double cluster::getUnWeightedX( const event::cluster_t &cluster, const Event& event){
  return event.dut->getPitchX() * cluster::getUnWeightedCol( cluster );
}

double cluster::getUnWeightedRow(const event::cluster_t &cluster){
  int n(0); double pos(0);
  for(int ii = 0; ii < cluster.size(); ii++){
    pos += cluster.at(ii)->row; n++;
  }
  return(pos/double(n));
}
double cluster::getUnWeightedY( const event::cluster_t &cluster, const Event& event){
  return event.dut->getPitchY() * cluster::getUnWeightedRow( cluster );
}

double cluster::getChargeWeightedCol(const event::cluster_t &cluster){
  double w(0), pos(0);
  for(int ii = 0; ii < cluster.size(); ii++){
    pos += cluster.at(ii)->col * cluster.at(ii)->tot;
    w += cluster.at(ii)->tot;
  }
  if(w == 0){
    cout << "cluster:: Cluster with 0tot detected. You have to check that this does not occur before calling getChargeWeightedCol" << endl;
    exit(-1);
  }
  return(pos/w);
}
double cluster::getChargeWeightedX( const event::cluster_t &cluster, const Event& event){
  return event.dut->getPitchX() * cluster::getChargeWeightedCol( cluster );
}

double cluster::getChargeWeightedRow(const event::cluster_t &cluster){
  double w(0), pos(0);
  for(int ii = 0; ii < cluster.size(); ii++){
    pos += cluster.at(ii)->row * cluster.at(ii)->tot;
    w += cluster.at(ii)->tot;
  }
  if(w == 0){
    cout << "cluster:: Cluster with 0tot detected. You have to check that this does not occur before calling getChargeWeightedCol" << endl;
    exit(-1);
  }
  return(pos/w);
}
double cluster::getChargeWeightedY( const event::cluster_t &cluster, const Event& event ){
  return event.dut->getPitchY() * cluster::getChargeWeightedRow( cluster );
}

double cluster::getEtaCorrectedCol(const event::cluster_t &cluster, DUT* dut){
  if (cluster.size() > 2){
    return -1;
  }else{
    double weighted = cluster::getChargeWeightedCol( cluster );
    double floored = floor(weighted);
    if(floored ==  weighted){ return weighted; }
    return( floored + dut->ecorrs.getX( weighted - floored) );
  }
}

double cluster::getEtaCorrectedX(const event::cluster_t &cluster, const Event& event, DUT* dut){
  return( event.dut->getPitchX() * cluster::getEtaCorrectedCol(cluster, dut) );
}


double cluster::getEtaCorrectedRow(const event::cluster_t &cluster, DUT* dut){
  if (cluster.size() > 2){
    return -1;
  }else{
    double weighted = cluster::getChargeWeightedRow( cluster );
    double floored = floor(weighted);
    if(floored ==  weighted){ return weighted; }
    return(floored + dut->ecorrs.getY( weighted - floored) );
  }
}


double cluster::getEtaCorrectedY(const event::cluster_t &cluster, const Event& event, DUT* dut){
  return( event.dut->getPitchY() * cluster::getEtaCorrectedRow(cluster, dut) );
}


double cluster::getMaxTotCellCol(const event::cluster_t &cluster){
  int maxCell = getMaxTotCell(cluster);
  return cluster.at(maxCell)->col;
}


double cluster::getMaxTotCellX(const event::cluster_t &cluster, const Event& event){
  return( event.dut->getPitchX() * getMaxTotCellCol(cluster) );
}


double cluster::getMaxTotCellRow(const event::cluster_t &cluster){
  int maxCell = getMaxTotCell(cluster);
  return cluster.at(maxCell)->row;
}


double cluster::getMaxTotCellY(const event::cluster_t &cluster, const Event& event){
  return( event.dut->getPitchY() * getMaxTotCellRow(cluster) );
}


bool cluster::isMatch(const event::cluster_t &cluster, const Event &event, int algo){
  bool match(false);
  int row = -1;
  int col = -1;
  
  if (algo == cluster::kUnweighted) {
    row = getUnWeightedRow(cluster);
    col = getUnWeightedCol(cluster);
  }else if (algo == cluster::kChargeWeighted){
    row = getChargeWeightedRow(cluster); 
    col = getChargeWeightedCol(cluster);
  }else if (algo == cluster::kEtaCorrected){
    row = getEtaCorrectedRow(cluster, event.dut); 
    col = getEtaCorrectedCol(cluster, event.dut);
  }else if (algo == cluster::kDigitalHeadTail1D){
    row = getClusterCenterRowByDigitalHeadTail1D(cluster); 
    col = getClusterCenterColByDigitalHeadTail1D(cluster);
  }else if (algo == cluster::kAnalogHeadTail1D){
    row = getClusterCenterRowByAnalogHeadTail1D(cluster, event); 
    col = getClusterCenterColByAnalogHeadTail1D(cluster, event);
  }else if (algo == cluster::kOneSidedAnalogHeadTail1D){
    row = getClusterCenterRowByOneSidedAnalogHeadTail1D(cluster, event);  
    col = getClusterCenterColByOneSidedAnalogHeadTail1D(cluster, event); 
  }
  
  if (algo!=cluster::kUnknown){
    if (row!=-1 && col!=-1){
      PllHit clusterPos;
      clusterPos.row = row;
      clusterPos.col = col;
      match = tbutils::isMatch( &clusterPos, event );
    }else{
      cerr << "cluster::isMatch: DUT:" << event.dut->getDUTid()
	  <<": Cluster position finding algorithm failed. Matching failed." << endl;
      match = false;
    }
  }else{
    for(int ii = 0; ii < cluster.size(); ii++){
      if( tbutils::isMatch( cluster.at(ii), event ) ){
	match = true;
	break;
      }
    }
  }
  
  return( match );
}

int cluster::getMatched(const vector< event::cluster_t > &clusters, const Event &event){
  int index(-1);
  for( int ii = 0; ii < clusters.size(); ii++){
    if( cluster::isMatch( clusters.at( ii ), event ) ){
      index = ii;
      break;
    }
  }
  return (index);
}

bool cluster::isCentralRegion(const event::cluster_t &cluster, const Event &event) {  
  for( int ii = 0; ii < cluster.size(); ii++) {
    const int col = cluster.at(ii)->col;
    const int row = cluster.at(ii)->row;
    if( col < event.dut->getSkipCols() ) return false;
    if( col > (event.dut->getMaxCol() - event.dut->getSkipCols()) ) return false;
    if( row < event.dut->getSkipRows() ) return false;
    if( row > (event.dut->getMaxRow() - event.dut->getSkipRows()) ) return false;
  }
  return true;
}

int cluster::spannedRows(const event::cluster_t &cluster){
  int min(999), max(-1);
  for(int ii =0; ii < cluster.size(); ii++){
    if( cluster.at(ii)->row < min) min = cluster.at(ii)->row;
    if( cluster.at(ii)->row > max) max = cluster.at(ii)->row;
  }
  return max - min + 1;
}

int cluster::spannedCols(const event::cluster_t &cluster){
  int min(999), max(-1);
  for(int ii =0; ii < cluster.size(); ii++){
    if( cluster.at(ii)->col < min) min = cluster.at(ii)->col;
    if( cluster.at(ii)->col > max) max = cluster.at(ii)->col;
  }
  return max - min + 1;
}

double cluster::getClusterCenterColByDigitalHeadTail1D(const event::cluster_t &cluster){
	int head,tail;
	//find head
	head = cluster.at(0)->col;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->col)<head){
			head = cluster.at(ii)->col;
		}
	}
	//find tail
	tail = cluster.at(0)->col;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->col)>tail){
			tail = cluster.at(ii)->col;
		}
	}
	return ((head+tail)/2);
}

double cluster::getClusterCenterXByDigitalHeadTail1D(const event::cluster_t &cluster, const Event& event ){
	  return( event.dut->getPitchX() * getClusterCenterColByDigitalHeadTail1D(cluster) );
}

double cluster::getClusterCenterRowByDigitalHeadTail1D(const event::cluster_t &cluster){
	int head,tail;
	//find head
	head = cluster.at(0)->row;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->row)<head){
			head = cluster.at(ii)->row;
		}
	}
	//find tail
	tail = cluster.at(0)->row;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->row)>tail){
			tail = cluster.at(ii)->row;
		}
	}
	return ((head+tail)/2);
}

double cluster::getClusterCenterYByDigitalHeadTail1D(const event::cluster_t &cluster, const Event& event ){
	  return( event.dut->getPitchY() * getClusterCenterRowByDigitalHeadTail1D(cluster) );
}

//FIXME: Lots of code doubling below

double cluster::getClusterCenterColByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	// signal_0: most probable signal released by a particle at 0 degree
	// t: 'detector thickness' or rather active thickness
	// theta: inclination
	// Beware that several presumptions (homogeneous charge collection in active zone,
	// active zone thickness, inclination angle...) are made that are not necessarily true.
  
	static long int warned = 0;
	if (!event.dut->getToTcalib()->hasCalib() || !event.dut->anglesCalculated){
		int dutBit = 1<<event.dut->getDUTid();
		if (!(warned & dutBit)){
		  cout << "cluster::getClusterCenterColByAnalogHeadTail1D: dut " << event.dut->getDUTid()
		    << ": missing calib or angles. Will not compute." << endl;
		  warned |= dutBit;
		}
		return -1;
	}
	
	// Check if angle is not too close to 90deg. We will divide by sin(eta) later
	double eta = event.dut->angleEtaCalculated;
	double etaTemp = abs(eta/pi);
	if (etaTemp - int(etaTemp) < 1E-6){
	  return -1;
	}
	
	double result = -1;
	double t = event.dut->getThickness();
	double signal_0 = event.dut->getSignal0();
	
	int head,tail;
	double signal_head_pixel, signal_tail_pixel;

	//find head
	//take into account when there are multiple head/tail-pixels
	head = cluster.at(0)->col;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->col)==head){
		  signal_head_pixel += event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}else if((cluster.at(ii)->col)<head){
			head = cluster.at(ii)->col;
			signal_head_pixel = event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}
	}
	
	//find tail
	tail = cluster.at(0)->col;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->col)==tail){
			signal_tail_pixel += event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}else if((cluster.at(ii)->col)>tail){
			tail = cluster.at(ii)->col;
			signal_tail_pixel = event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);	
		}
	}

	double signal_theta;
	signal_theta = (signal_0 * event.dut->getPitchX())/(t * sin(eta));

	result = getClusterCenterColByDigitalHeadTail1D(cluster)
	      +(event.dut->getPitchX())*(min((double) signal_head_pixel,signal_theta)
	      -min((double) signal_tail_pixel,signal_theta))/(2*signal_theta);

	return result;
}

double cluster::getClusterCenterXByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	return( event.dut->getPitchX() * getClusterCenterColByAnalogHeadTail1D(cluster, event) );
}

double cluster::getClusterCenterRowByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	// signal_0: most probable signal released by a particle at 0 degree
	// t: 'detector thickness' or rather active thickness
	// phi: inclination
	// Beware that several presumptions (homogeneous charge collection in active zone,
	// active zone thickness, inclination angle...) are made that are not necessarily true.
	
	static long int warned = 0;
	if (!event.dut->getToTcalib()->hasCalib() || !event.dut->anglesCalculated){
		int dutBit = 1<<event.dut->getDUTid();
		if (!(warned & dutBit)){
		  cout << "cluster::getClusterCenterRowByAnalogHeadTail1D: dut " << event.dut->getDUTid()
		    << ": missing calib or angles. Will not compute." << endl;
		  warned |= dutBit;
		}
		return -1;
	}
	
	// Check if angle is not too close to 90deg. We will divide by sin(phi) later
	double phi = event.dut->anglePhiCalculated;
	double phiTemp = abs(phi/pi);
	if (phiTemp - int(phiTemp) < 1E-6){
	  return -1;
	}
	
	double result = -1;
	double t = event.dut->getThickness();
	double signal_0 = event.dut->getSignal0();

	int head = -1,tail = -1;
	double signal_head_pixel = 0, signal_tail_pixel = 0;

	//find head
	head = cluster.at(0)->row;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->row)==head){
		  signal_head_pixel += event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}else if((cluster.at(ii)->row)<head){
			head = cluster.at(ii)->row;
			signal_head_pixel = event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, head);
		}
	}
	//find tail
	tail = cluster.at(0)->row;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->row)==tail){
			signal_tail_pixel += event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}else if((cluster.at(ii)->row)>tail){
			tail = cluster.at(ii)->row;
			signal_tail_pixel = event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, tail);	
		}
	}

	double signal_phi;
	signal_phi = (signal_0 * event.dut->getPitchY())/(t * sin(phi));

	result = getClusterCenterRowByDigitalHeadTail1D(cluster)+(event.dut->getPitchY())
	      *(min((double) signal_head_pixel,signal_phi)
		-min((double) signal_tail_pixel,signal_phi))
	      /(2*signal_phi);
	return result;
}

double cluster::getClusterCenterYByAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	return( event.dut->getPitchY() * getClusterCenterRowByAnalogHeadTail1D(cluster, event) );
}

double cluster::getClusterCenterColByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	// signal_0: most probable signal released by a particle at 0 degree
	// t: 'detector thickness' or rather active thickness
	// theta: inclination
	// Beware that several presumptions (homogeneous charge collection in active zone, 
        // active zone thickness, inclination angle...) are made that are not necessarily true.
	
	static long int warned = 0;
	if (!event.dut->getToTcalib()->hasCalib() || !event.dut->anglesCalculated){
		int dutBit = 1<<event.dut->getDUTid();
		if (!(warned & dutBit)){
		  cout << "cluster::getClusterCenterColByOneSidedAnalogHeadTail1D: dut " << event.dut->getDUTid()
		    << ": missing calib or angles. Will not compute." << endl;
		  warned |= dutBit;
		}
		return -1;
	}
	
	// Check if angle is not too close to 90deg. We will divide by sin(eta) later
	double eta = event.dut->angleEtaCalculated;
	double etaTemp = abs(eta/pi);
	if (etaTemp - int(etaTemp) < 1E-6){
	  return -1;
	}
	
	double t = event.dut->getThickness();
	double signal_0 = event.dut->getSignal0();  
	
	double result = -1;
	int head = -1;
	int signal_head_pixel = 0;

	//find head
	// FIXME: Determine which side is supposed to be the head side.
	head = cluster.at(0)->col;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->col)==head){
		  signal_head_pixel += event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}else if((cluster.at(ii)->col)<head){
			head = cluster.at(ii)->row;
			signal_head_pixel = event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}
	}
	

	double signal_theta;
	signal_theta = (signal_0 * event.dut->getPitchX())/(t * sin(eta));

	result = head + 0.5 * event.dut->getPitchX() - event.dut->getPitchX() 
	    * signal_head_pixel/signal_theta + 0.5 * t * tan(eta);

	return result;
}

double cluster::getClusterCenterXByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	return( event.dut->getPitchX() * getClusterCenterColByOneSidedAnalogHeadTail1D(cluster, event) );
}

double cluster::getClusterCenterRowByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	// signal_0: most probable signal released by a particle at 0 degree
	// t: 'detector thickness' or rather active thickness
	// theta: inclination
	// Beware that several presumptions (homogeneous charge collection in active zone, 
	// active zone thickness, inclination angle...) are made that are not necessarily true.
	
	static long int warned = 0;
	if (!event.dut->getToTcalib()->hasCalib() || !event.dut->anglesCalculated){
		int dutBit = 1<<event.dut->getDUTid();
		if (!(warned & dutBit)){
		  cout << "cluster::getClusterCenterRowByOneSidedAnalogHeadTail1D: dut " << event.dut->getDUTid()
		    << ": missing calib or angles. Will not compute." << endl;
		  warned |= dutBit;
		}
		return -1;
	}
	
	// Check if angle is not too close to 90deg. We will divide by sin(phi) later
	double phi = event.dut->anglePhiCalculated;
	double phiTemp = abs(phi/pi);
	if (phiTemp - int(phiTemp) < 1E-6){
	  return -1;
	}

	double t = event.dut->getThickness();
	double signal_0 = event.dut->getSignal0();  
	
	double result = -1;
	int head = -1;
	int signal_head_pixel = 0;

	//find head
	// FIXME: Determine which side is supposed to be the head side.
	head = cluster.at(0)->row;
	for(int ii = 0; ii < cluster.size(); ii++){
		if((cluster.at(ii)->row)==head){
		  signal_head_pixel += event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}else if((cluster.at(ii)->row)<head){
			head = cluster.at(ii)->row;
			signal_head_pixel = event.dut->q(cluster.at(ii)->tot,
				      cluster.at(ii)->col, cluster.at(ii)->row);
		}
	}

	double signal_phi;
	signal_phi = (signal_0 * event.dut->getPitchY())/(t * sin(phi));

	result = head + 0.5 * event.dut->getPitchY() - event.dut->getPitchY() 
	      * signal_head_pixel/signal_phi + 0.5 * t * tan(phi);

	return result;
}


double cluster::getClusterCenterYByOneSidedAnalogHeadTail1D(const event::cluster_t &cluster, const Event& event){
	return( event.dut->getPitchY() * getClusterCenterRowByOneSidedAnalogHeadTail1D(cluster, event) );
}


double cluster::getCol(const event::cluster_t &cluster, const Event& event, int algo)
{
  if (algo == kUnweighted) {
    return getUnWeightedCol(cluster);
  }else if (algo == kChargeWeighted){
    return getChargeWeightedCol(cluster);
  }else if (algo == kEtaCorrected){
    return getEtaCorrectedCol(cluster, event.dut);
  }else if (algo == kMaxCell){
    return getMaxTotCellCol(cluster);
  }else if (algo == kDigitalHeadTail1D){
    return getClusterCenterColByDigitalHeadTail1D(cluster);
  }else if (algo == kAnalogHeadTail1D){
    return getClusterCenterColByAnalogHeadTail1D(cluster, event);
  }else if (algo == kOneSidedAnalogHeadTail1D){
    return getClusterCenterColByOneSidedAnalogHeadTail1D(cluster, event); 
  }else{
    static long int warned = 0;
    if (! (warned & 1<<algo)){
      cerr << "cluster::GetCol: Algorithm " << algo << " unknown" << endl;
      warned |= 1<<algo;
    }
    return -1;
  }
}


double cluster::getRow(const event::cluster_t &cluster, const Event& event, int algo)
{
  if (algo == kUnweighted) {
    return getUnWeightedRow(cluster);
  }else if (algo == kChargeWeighted){
    return getChargeWeightedRow(cluster);
  }else if (algo == kEtaCorrected){
    return getEtaCorrectedRow(cluster, event.dut);
  }else if (algo == kMaxCell){
    return getMaxTotCellRow(cluster);
  }else if (algo == kDigitalHeadTail1D){
    return getClusterCenterRowByDigitalHeadTail1D(cluster);
  }else if (algo == kAnalogHeadTail1D){
    return getClusterCenterRowByAnalogHeadTail1D(cluster, event);
  }else if (algo == kOneSidedAnalogHeadTail1D){
    return getClusterCenterRowByOneSidedAnalogHeadTail1D(cluster, event);
  }else{
    static long int warned = 0;
    if (! (warned & 1<<algo)){
      cerr << "cluster::GetRow: Algorithm " << algo << " unknown" << endl;
      warned |= 1<<algo;
    }
    return -1;
  }
}


double cluster::getX(const event::cluster_t &cluster, const Event& event, int algo)
{
  if (algo == kUnweighted) {
    return getUnWeightedX( cluster, event );
  }else if (algo == kChargeWeighted){
    return getChargeWeightedX( cluster, event );
  }else if (algo == kEtaCorrected){
    return getEtaCorrectedX(cluster, event, event.dut);
  }else if (algo == kMaxCell){
    return getMaxTotCellX(cluster, event);
  }else if (algo == kDigitalHeadTail1D){
    return getClusterCenterXByDigitalHeadTail1D(cluster, event );    
  }else if (algo == kAnalogHeadTail1D){
    return getClusterCenterXByAnalogHeadTail1D(cluster, event);
  }else if (algo == kOneSidedAnalogHeadTail1D){
    return getClusterCenterXByOneSidedAnalogHeadTail1D(cluster, event);
  }else{
    static long int warned = 0;
    if (! (warned & 1<<algo)){
      cerr << "cluster::GetCol: Algorithm " << algo << " unknown" << endl;
      warned |= 1<<algo;
    }
    return -1;
  }
}


double cluster::getY(const event::cluster_t &cluster, const Event& event, int algo)
{
  if (algo == kUnweighted) {
     return getUnWeightedY( cluster, event );
  }else if (algo == kChargeWeighted){
     return getUnWeightedY( cluster, event );
  }else if (algo == kEtaCorrected){
    return getEtaCorrectedY(cluster, event, event.dut);
  }else if (algo == kMaxCell){
    return getMaxTotCellY(cluster, event);
  }else if (algo == kDigitalHeadTail1D){
    return getClusterCenterYByDigitalHeadTail1D(cluster, event );
  }else if (algo == kAnalogHeadTail1D){
    return getClusterCenterYByAnalogHeadTail1D(cluster, event);
  }else if (algo == kOneSidedAnalogHeadTail1D){
    return getClusterCenterYByOneSidedAnalogHeadTail1D(cluster, event);
  }else{
    static long int warned = 0;
    if (! (warned & 1<<algo)){
      cerr << "cluster::GetCol: Algorithm " << algo << " unknown" << endl;
      warned |= 1<<algo;
    }
    return -1; 
  }
}
  



