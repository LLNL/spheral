// managed.cpp
#include <iostream>

#include "../Timer.hh"

Timer Everything("Everything");
Timer loop1(" Loop 1", Everything);
Timer loop2(" Loop 2", Everything);

int main() {

  Everything.start();
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
  Everything.stop();

  Timer::TimerSummary();
}

