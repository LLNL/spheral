// unmanaged.cpp
#include <iostream>

#include "../Timer.hh"

Timer loop1, loop2;

int main() {

  double y1=0, y2=0;

  loop1.start();
  for(int i=0; i<12345678; i++) {
    y1 += (float)i;
  }
  loop1.stop();

  loop2.start();
  for(int i=0; i<12345678; i++) {
    y2 += (float)i/(i+1);
  }
  loop2.stop();

  cout << " y1= " << y1 << " y2= " << y2 << endl;
  cout << " loop1 time=" << loop1.wc_time() 
       << " loop2 time=" << loop2.wc_time() << endl;

}

