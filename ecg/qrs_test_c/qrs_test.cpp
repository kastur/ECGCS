#include <iostream>
#include "/home/k/ecg/third/eplimited/qrsdet2-inl.h"
#include "/home/k/ecg/third/eplimited/qrsfilt-inl.h"
#include "/home/k/ecg/mbed/analog-in-inl.h"

using namespace std;
using qrsdet::QRSDet;
using qrsdet::deriv1;

void detect_run(const string&);

int main(int argc, char** argv) {
  if (argc != 2) {
    cerr << "Provide ecg data file as argument" << endl;
    return -1;  
  }
  string file = argv[1];
  detect_run(file);
  return 0;
}

void detect_run(const string& file) {
  AnalogIn ecg_in(file, 1.0, 0.0);
  for (int i = 0; i < 5000; i++) {
    float float_val = ecg_in.read();
    int int_val = float_val * 200;
    int init = (i == 0) ? 1 : 0;
    int computed_val = qrsdet::QRSDet(int_val, init);
    cout << int_val << "\t" << computed_val << endl;
  }
}
