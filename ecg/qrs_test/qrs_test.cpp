#include <iostream>
#include "/home/k/ecg/third/eplimited/qrsdet2-inl.h"
#include "/home/k/ecg/third/eplimited/qrsfilt-inl.h"
#include "/home/k/ecg/mbed/analog-in-inl.h"
#include "/home/k/ecg/usrlib/usrlib.h"

using namespace std;
using qrsdet::QRSDet2;

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
  QRSDet2<8> qrs_det;
  qrsdet::DelayFilter<256> ecg_buffer;
  qrsdet::DelayFilter<256> qrs_buffer;
  for (int i = 0; i < 5000; i++) {
    float float_val = ecg_in.read();
    int int_val = float_val * 200;
    int delay_ecg = ecg_buffer.Process(int_val);
    int delay_qrs = qrs_buffer.Process(0);
    cout << delay_ecg << "\t" << delay_qrs << "\n";
    int qrs_val= qrs_det.Process(int_val);
    if (qrs_val > 0) {
      for (int i = 0; i < 10; i++)
        qrs_buffer.SetValue(100, (-qrs_val+i));
/*
      int n;
      int* morphology = 
        buffer.GetRange(-(qrs_val+ 50), -(qrs_val - 45), &n);
      for(int i = 0; i < n; i++) {
        cout << morphology[i] << "\t";
      }
      cout << "\n";
      delete [] morphology;
*/
    }
  }
}
