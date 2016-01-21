// Timer.hh
#ifndef TIMER_H
#define TIMER_H

#ifdef USE_MPI
#define MPI
#endif

#ifdef MPI
#include "mpi.h"
#include "Distributed/Communicator.hh"
#endif

// if TIMER is not defined, then we just have a stub class.
#ifdef TIMER

extern "C" {
  // these are for gettimeofday() wall-clock timing
#include <unistd.h>
#include <sys/time.h>

  // obviously, if you don't have PAPI installed on your system
  // (which requires kernel modification) then the counters wont 
  // work, so don't compile with -DPAPI
#ifdef PAPI
#include <papi.h>
#endif

}

#ifndef TIMER_COUNTER
#define TIMER_COUNTER 0
#endif

#include <list>
#include <string>
using namespace std;

#define DIAGNOSTIC false 

class Timer {
  
public:
  
  Timer();                             // unmanaged
  Timer(const string&);                // root
  Timer(const string&, Timer&);        // managed
  Timer(const string&, Timer&, bool);  // diagnostic
  ~Timer();
  
  void setup();
  void start();
  void stop();
  void clear();

  inline double getTimeStampWC();
  double wc_time() { return accumulated_WCtime; }

#ifdef PAPI
  long long int papi_counter1() { return accumulated_counter1; }
  long long int papi_counter2() { return accumulated_counter2; }
#endif

  string Name() { return timer_name;  }
  
  long int Count() { return count; }

  static list<Timer*> TimerList;

  static void TimerSummary(const int bert, const int ernie) {
    TimerSummary(); // backwards compatibilty...
  }
  
  static void TimerSummary(void);
  
private:
  
  //bool timer_on;    // State of timer, either on(true) or off(false)
  bool diagnostic;

  double accumulated_WCtime, last_WCtime_stamp;

  // wall-clock timer data
#ifndef MPI  
  struct timeval tv;   //  Values from call to gettimeofday
  struct timezone tz;
#endif

  int ID;  
  string timer_name;
  Timer& Parent;

  double avgWC,  minWC,  maxWC;

#ifdef PAPI
  long long int values[2];
  long long int accumulated_counter1;
  long long int accumulated_counter2;

#endif

  long long int count;
};


#else


// stub Timer class
#include <string>
#include <iostream>
using namespace std;

#define DIAGNOSTIC false

class Timer {

public:

  Timer() {}
  Timer(const string&) {}
  //Timer(const string&, int) {}
  Timer(const string&, Timer&) {}
  Timer(const string&, Timer&, bool) {}
  ~Timer() {}

  //static list<Timer*> TimerList;

  inline void setup(){}
  inline void start(){}
  inline void stop(){}  
  inline void clear(){}

  inline double getTimeStampWC() {return 0.0;}
  inline double wc_time() {return 0.0;}

#ifdef PAPI
  inline long long int papi_counter1() {return 0;}
  inline long long int papi_counter2() {return 0;}
#endif

  inline string Name() {return NULL;}
  
  inline long int Count() {return 0;}
  
  static void TimerSummary(const int bert, const int ernie) {
    TimerSummary(); // backwards compatibilty...
  }

  static void TimerSummary(void) {
    int rank;
#ifdef MPI
    MPI_Comm_rank(Communicator::communicator(), &rank);
#else
    rank=0;
#endif
    if(rank==0) {
      cout << " Timers Disabled.  No timing output written."  << endl;
    }
  }

};

#endif // TIMER

#endif // TIMER_H
