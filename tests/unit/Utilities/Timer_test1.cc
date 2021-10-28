// Code that tests the Timer classes

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "../Timer.hh"
#include "Distributed/Communicator.hh"

#define N0 10
#define N1 90
#define N2 900
#define N3 9000
#define M  1000.0

// Must initialize the static list defined in Timer.hh
// i don't know how to get rid of this...
list<Timer*> Timer::TimerList(0); 

// Root Timer declaration.  A Timer that is declared without passing a
// parent Timer is used as the root of the tree of Timers.  Currently,
// there is no check to ensure that the user only declares one root
// Timer...  The Everything Timer times the entire code.
Timer Everything("Everything");

// By passing a reference of Everything Timer to the contructor, the
// TimerInfo class will know that the Work Timer is part of or inside
// of the parent (Everything)
Timer Work      ("Work",     Everything);

// The above Timers are defined globally because there are Timers
// in the goofy_function that need access to the Work Timer.
void goofy_function(int);


// Diagnostic Timers.  When it's not meaningful to specify a
// parent Timer (or you can't -- or don't care) you can create
// a Diagnostic Timer that behaves as normal Timers, but will
// be printed out in a separate table in the output.  The
// parent specified in the construction is only to compute
// table percentages.  The 3rd arguement in the contruction 
// is what makes this a "diagnostic" Timer.
Timer Diag1("Diag1", Everything, DIAGNOSTIC);
Timer Diag2("Diag2", Everything, DIAGNOSTIC);

int main(int argc, char **argv)  {
  
  int rank, number_procs;
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(Spheral::Communicator::communicator(), &rank);
  MPI_Comm_size(Spheral::Communicator::communicator(), &number_procs);
#else
  rank=0;
  number_procs=1;
#endif
  cout << "rank = " << rank << " number_procs= " << number_procs << endl;


  //  An Unmanaged Timer is declared as follows:
  Timer Unmanaged;

  // The Unmanaged Timer can be used as a regular Timer, but will
  // not go into the Summary.  Therefore the user must make the 
  // appropriate calls to obtain the accumulated time.
  Unmanaged.start();
  // do something
  Unmanaged.stop();
  cout << " Unmanaged Time:  " << Unmanaged.wc_time() << endl;

  Timer o;
  
  o.start();
  Timer Create    ("Create",   Everything);
  Timer Stuff     ("Stuff",    Everything);

  Timer RAM    ("RAM",     Work);
  Timer Barrier("Barrier", Work);
  o.stop();
  cout << " o wc_time=" << o.wc_time() 
       << " count =" << o.Count() << endl;

  o.start();
  o.stop();
  cout << " o wc_time=" << o.wc_time() 
       << " count =" << o.Count() << endl;


  Everything.start();

  Create.start();
  int N = (int)5e5;
  vector<double> big;

  cout << "N = " << N << endl;

  for(int i=0; i< N ; ++i) {

    big.push_back(i*1.0);

  }
  Create.stop();

  cout << " rank" << rank << "   Create wc= " << Create.wc_time() << endl;  

  double r1, sum;
  srand(rank);
  
  Work.start();
  
  for(int j=1 ; j<=6 ; j++) {

    sum = 0.0;

    RAM.start();
    for (int i=0; i<N/10 ; ++i) {
      r1 = (double(rand())/RAND_MAX)*N;
      
      sum += big[(int)r1];
      
      //cout << "i: " << i << " r1: " << r1 << " big[r1]: " << big[r1] << endl;
    }
    RAM.stop();

    Barrier.start();
#ifdef USE_MPI
    MPI_Barrier(Spheral::Communicator::communicator());
#endif
    Barrier.stop();

  }

  goofy_function(N3);

  Diag1.start();
  goofy_function(N2);
  Diag1.stop();

  Diag2.start();
  goofy_function(N1);
  Diag2.stop();

  Work.stop();

  Everything.stop();

  Timer::TimerSummary();
  
#ifdef USE_MPI
  MPI_Finalize();
#endif
}


void goofy_function(int N) {
  
  static  Timer Chew("Chew", Work);
  
  static Timer A1("A1", Chew);
  static Timer A2("A2", Chew);
  static Timer A3("A3", Chew);
  static Timer A4("A4", Chew);
  
  Chew.start();
  
  double y;
  double k=0;
  for ( double x=1.; x<=N1; x+=1) {
    
    A1.start();
    for(y=1.0; y<=M; y+=.1) k += y;
    A1.stop();
    
    A2.start();
    for(y=1.0; y<=M; y+=.1) k += sin(y)/y;
    A2.stop();
    
    A3.start();
    for(y=1.0; y<=M; y+=.1) k += sin(y) + y;
    A3.stop();
    
    A4.start();
    A4.stop();
  }
  
  cout << " k = " << k << endl;
  
  Chew.stop();
 
}

