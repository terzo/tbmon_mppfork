#include "simThreeVector.h"
#include "event.h"

//Partly copypasta from
// http://bytes.com/topic/c/insights/652951-how-parse-file-c

istream& operator>>(istream& is, simThreeVector& vec) throw (ios_base::failure) {
  string tmp;

  ios_base::iostate oldIOState = is.exceptions();
  cin.exceptions(~ios::goodbit);  // turn on exceptions
  try {
    getline (is, tmp, '('); //Throw prepending stuff
    getline (is, tmp, ',');
    vec.data[0] = atof(tmp.c_str());
    getline (is, tmp, ',');
    vec.data[1] = atof(tmp.c_str());
    getline (is, tmp, ')');
    vec.data[2] = atof(tmp.c_str());
  }
  catch(ios_base::failure failure) {
    assert(!is.good());
    cout << "[ simThreeVector ]; Some error in parsing!" << endl;
    exit(-1);
    is.exceptions(oldIOState);  // restoring old IO exception handling
    throw; // there is no way to recover the stream without more info
  }
  is.exceptions(oldIOState);  // restoring old IO exception handling
  return is;
}

string simThreeVector::toString() {
  ostringstream out;
  out << "(" << data[0] << "," << data[1] << "," << data[2] << ")";
  return out.str();
}

//Translate from sensor-local coordinate system
// as used in TestBeamSim -> what is used in tbmon
void simThreeVector::localTranslate(const Event& event) {
  double tmp = data[0];
  //data[0] = data[1]*1000 + (event::ncols/2.0 - 0.5)*event::pitchX;
  //data[1] = tmp*1000 + (event::nrows/2.0 - 0.5)*event::pitchY;
  data[0] = data[1]*1000 + (event.dut->getNcols()/2.0 - 0.5)*event.dut->getPitchX();
  data[1] = tmp*1000 + (event.dut->getNrows()/2.0 - 0.5)*event.dut->getPitchY();
}
