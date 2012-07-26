#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void periodic_solver(Group& group, Fractal_Memory& mem,Fractal& frac)
  {
    cout << "entering periodic " << endl;
    dens_to_slices(group,mem,frac);
    const int length_1=frac.get_grid_length();  
    const int length_c=(length_1+2)/2;
    const int length_2=length_1*length_1;
    const double length_3=(double)(length_1*length_2);
    const double pi=4.0*atan(1.0);
    group.set_force_const(4.0*pi/static_cast<double>(length_2));
    size_t sizeR=sizeof(double);
    size_t sizeC=sizeof(fftw_complex);
    mem.p_mess->potR=(double*) fftw_malloc(2*sizeR*mem.p_mess->total_memory);
    mem.p_mess->potC=(fftw_complex*) fftw_malloc(sizeC*mem.p_mess->total_memory);
    vector <double> green_1(length_1+1);
    for (int k=0;k <length_1+1;++k)
      {
	double aa=pi*(double)(k)/(double)(length_1);
	green_1[k]=1.0e-30+pow(2.0*sin(aa),2);
      }
    receive_dens(mem,frac,length_1);
    fftw_execute(mem.p_mess->plan_rc);
    double g_c=group.get_force_const()/length_3;
    cout << "going to power spectrum from periodic " << length_1 << endl;
    vector <double> variance_rho(length_1,0.0);
    vector <double> variance_pot(length_1,0.0);
    vector <double> variance_force(length_1,0.0);
    vector <double> variance_force_s(length_1,0.0);
    if(mem.p_mess->start_x == 0)
      {
	mem.p_mess->potC[0][0]=0.0;
	mem.p_mess->potC[0][1]=0.0;
      }
    frac.timing(-1,6);
    power_spectrum(mem.p_mess->potC,length_1,variance_rho,variance_pot,variance_force,variance_force_s,0,frac.get_density_0(),mem.do_var,mem);

    frac.timing(1,6);
    for(int nx=0;nx < mem.p_mess->length_x;++nx)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nz=0;nz < length_c;++nz)
	      {
		int twit=mem.fftw_where(nx,ny,nz,length_1,length_c);
		double of_the_year=-g_c/(green_1[nz]+green_1[ny]+green_1[nx+mem.p_mess->start_x]);
		mem.p_mess->potC[twit][0]*=of_the_year;
		mem.p_mess->potC[twit][1]*=of_the_year;
	      }
	  }
      }
    fftw_execute(mem.p_mess->plan_cr);
    slices_to_potf(mem,frac,length_1);
    fftw_free(mem.p_mess->potR);
    fftw_free(mem.p_mess->potC);
    receive_potf(group,mem,frac);
    cout << "exiting periodic " << endl;
  }
}



