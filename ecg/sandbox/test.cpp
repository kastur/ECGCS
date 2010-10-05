#include <iostream>
using namespace std;


unsigned int random_number() {
	static unsigned int lfsr = 1;
	static unsigned int period = 0;
	do {
	  /* taps: 32 31 29 1; characteristic polynomial: x^32 + x^31 + x^29 + x + 1 */
	  lfsr = (lfsr >> 1) ^ (unsigned int)(0 - (lfsr & 1u) & 0xd0000001u); 
	  ++period;
	} while(lfsr != 1u);
	return lfsr;
}

int main() {
  cout << random_number() << endl;
  cout << random_number() << endl;
  cout << random_number() << endl;
  cout << random_number() << endl;

  return 0;
}
