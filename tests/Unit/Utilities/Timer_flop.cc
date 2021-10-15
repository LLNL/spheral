// flop.cc
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "../Timer.hh"

#define N0 1
#define NLOOP 1000000


// Must initialize the static list defined in Timer.hh
//list<Timer*> Timer::TimerList(0); 

void flopping_test(void);


Timer Everything("Everything");
Timer TotalFlopping("TotalFlopping",     Everything);
Timer flopalloc ("flop alloc",TotalFlopping);
Timer flop0     ("flop0",    TotalFlopping);
Timer flop1     ("flop1",    TotalFlopping);
Timer flop2     ("flop2",    TotalFlopping);
Timer flop3     ("flop3",    TotalFlopping);
Timer flop4     ("flop4",    TotalFlopping);
Timer flop5     ("flop5",    TotalFlopping);
Timer flop6     ("flop6",    TotalFlopping);
Timer flop7     ("flop7",    TotalFlopping);
Timer flop8     ("flop8",    TotalFlopping);
Timer flop9     ("flop9",    TotalFlopping);

int main() {
  Everything.start();
  int rank, number_procs;
  rank=0;
  number_procs=1;
  
  TotalFlopping.start();
  for(int i=0; i<N0; i++) {
    flopping_test();
  }
  TotalFlopping.stop();


  Everything.stop();
  Timer::TimerSummary();  
}

void flopping_test(void) {

  flopalloc.start();
  int i;
  double y=0, f=0;
  double a[NLOOP];
  flopalloc.stop();

  flop0.start();  
  for(i=0; i<NLOOP; i++) a[i] = 0.0;
  flop0.stop();

  flop1.start();  
  for(i=0; i<NLOOP; i++) a[i] = 1.01;
  flop1.stop();

  flop2.start();
  for(i=0; i<NLOOP; i++) y+=i;  // 2 FP
  flop2.stop();

  flop3.start();
  for(i=0; i<NLOOP; i++) y+=3.2;  // 1 FP
  flop3.stop();

  flop4.start();
  for(i=0; i<NLOOP; i++) a[i] = y + i;           // 2 FP
  flop4.stop();

  flop5.start();
  for(i=0; i<NLOOP; i++) a[i] = a[i] * y;       // 1 FP
  flop5.stop();

  flop6.start();
  for(i=0; i<NLOOP; i++) y += a[i];              // 1 FP
  flop6.stop();

  flop7.start();
  for(i=0; i<NLOOP; i++) a[i] = a[i] * y + a[i];  // 2 FP
  flop7.stop();


  flop8.start();
  for(i=0; i<NLOOP; i++) {f=i*y; y+=f; f=y*f;}
  flop8.stop();

  flop9.start();
  for(i=0; i<NLOOP; i++) f=f/(i+1);
  flop9.stop();



  printf(" y=%e  f=%e  %e  %e\n", y, f, a[0], a[NLOOP-1]);

}
