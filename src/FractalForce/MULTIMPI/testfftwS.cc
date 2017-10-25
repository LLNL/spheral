#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <algorithm>
#include <complex>
#include <ctime>
#include <climits>
#include <cerrno>
#include <sys/stat.h>
using namespace std;
#include <fftw3.h>
int main(int argc, char **argv){
  const ptrdiff_t N0 = 32, N1 = 32;
  
  fftw_plan plan1,plan2;
  double *data1in; //local data of course
  fftw_complex *data1out; //local data of course
  double *data2in; //local data of course
  fftw_complex *data2out; //local data of course
  ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
  
  /* get local data size and allocate */
  alloc_local = N0*(N1/2+1);
  data1out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
  data1in  = (double *) fftw_malloc(2*sizeof(double) * alloc_local);
  data2out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
  data2in  = (double *) fftw_malloc(2*sizeof(double) * alloc_local);
  /* create plan for forward DFT */
  plan1 = fftw_plan_dft_r2c_2d(N0, N1, data1in, data1out,FFTW_ESTIMATE);
  plan2 = fftw_plan_dft_r2c_2d(N0, N1, data2in, data2out,FFTW_ESTIMATE);
  
  /* initialize data to some function my_function(x,y) */
  int N1A=(N!/2+1);
  int N1B=2*N1A;
  for (i = 0; i < local_n0; ++i) 
    {
      for (j = 0; j < N1; ++j)
	{
	  data1in[i*N1B+j]=1.0;
	  data2in[i*N1B+j]=0.0;
	}
    }
  data2in[0]=1.0;
  /* compute transforms, in-place, as many times as desired */
  fftw_execute(plan1);
  fftw_execute(plan2);
  for (i = 0; i < N0; ++i) for (j = 0; j < N1A; ++j){
    cout << " out " << local_0_start << " " << i << " " << j << " " << data1out[i*N1A+j][0] << " " << data1out[i*N1A+j][1];
    cout << " " << data2out[i*N1A+j][0] << " " << data2out[i*N1A+j][1] << endl;
  }
  
  fftw_destroy_plan(plan1);
  fftw_free(data1out);
  fftw_free(data1in);
  fftw_destroy_plan(plan2);
  fftw_free(data2out);
  fftw_free(data2in);
  printf("finalize\n");
  return 0;
}
