// diagnostic.cpp
#include <iostream>

#include "../Timer.hh"

Timer Everything("Everything");
Timer loop1(" Loop 1", Everything);
Timer loop2(" Loop 2", Everything);

double do_work1(void);
double do_work2(void);
double transform(double y);

int main() {

  Everything.start();
  double y1=0, y2=0;

  loop1.start();
  y1 = do_work1();
  loop1.stop();

  loop2.start();
  y2 = do_work2();
  loop2.stop();

  cout << " y1= " << y1 << " y2= " << y2 << endl;
  Everything.stop();

  Timer::TimerSummary();
}


double do_work1() {
  double y1;
  for(int i=0; i<1234567; i++) {
    y1 += (double)i;
  }

  y1 = transform(y1);
  return y1;
}

double do_work2() {
  double y2;
  for(int i=0; i<12345678; i++) {
    y2 += (double)i/(i+1.0);
  }
  y2 = transform(y2);
  return y2;
}


double transform(double y) {
  static Timer Transform("Transform", Everything, DIAGNOSTIC);

  Transform.start();
  int N = (int)(y/12345678);

  cout << y << "   " << N << endl;

  for(int i=0; i<N; i++) {
    y -= (double)i;
  }
  Transform.stop();

  return y;
}
