#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    cout << "entering isolated " << endl;
    dens_to_slices(group,mem,frac);
    const int length_1=frac.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    ptrdiff_t total_memory,start_x,length_x;
    if(Fractal::first_time_solver)
      {
	Fractal::first_time_solver=false;
	return;
      }
    //
    cout << "isol 0 " << endl;
    total_memory=fftw_mpi_local_size_3d(length_2,length_2,length_11,MPI::COMM_WORLD,&length_x,&start_x);
    size_t sizeR=sizeof(double);
    size_t sizeC=sizeof(fftw_complex);
    mem.p_mess->potR=(double*) fftw_malloc(2*sizeR*mem.p_mess->total_memory);
    mem.p_mess->potC=(fftw_complex*) fftw_malloc(sizeC*mem.p_mess->total_memory);
    mem.p_mess->zeroR();
    cout << "isol a " << endl;
    receive_dens(mem,frac,length_2);
    fftw_execute(mem.p_mess->plan_rc);
    //
    for(int n_x=0;n_x < length_x;++n_x)
      {
	int n_a=min(n_x,length_2-n_x);
	for(int n_y=0;n_y < length_2;++n_y)
	  {
	    int n_b=min(n_y,length_2-n_y);
	    for(int n_z=0;n_z < length_11;++n_z)
	      {
		int n_c=min(n_z,length_2-n_z);
		int wherexyz=mem.fftw_where(n_x,n_y,n_z,length_2,length_11);
		int whereabc=mem.fftw_where(n_a,n_b,n_c,length_11,length_11);
		mem.p_mess->potC[wherexyz][0]*=mem.p_mess->green[whereabc];
		mem.p_mess->potC[wherexyz][1]*=mem.p_mess->green[whereabc];
	      }
	  }
      }
    cout << "isol b " << endl;
    fftw_execute(mem.p_mess->plan_cr);
    cout << "isol c " << endl;
    slices_to_potf(mem,frac,length_2);
    cout << "isol d " << endl;
    fftw_free(mem.p_mess->potR);
    fftw_free(mem.p_mess->potC);
    receive_potf(group,mem,frac);
    cout << "exiting isolated " << endl;
  }
}
