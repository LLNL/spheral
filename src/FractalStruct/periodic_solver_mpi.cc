#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void periodic_solver(Group& group, Fractal_Memory& mem,Fractal& frac)
  {
    frac.timing(-1,4);
    ofstream& FileFFT=frac.p_file->DUMPS;
    FileFFT << "entering periodic " << "\n";
    const int length_1=frac.get_grid_length();  
    const int length_c=(length_1+2)/2;
    const int length_2=length_1*length_1;
    const double length_3=static_cast<double>(length_1*length_2);
    const double pi=4.0*atan(1.0);
    group.set_force_const(4.0*pi/static_cast<double>(length_2));
    vector <double> green_1(length_1+1);
    for (int k=0;k <length_1+1;++k)
      {
	double aa=pi*static_cast<double>(k)/static_cast<double>(length_1);
	green_1[k]=1.0e-30+pow(2.0*sin(aa),2);
      }
    Full_Stop(mem,34);
    frac.timing(-1,24);
    dens_to_slices(group,mem,frac);
    frac.timing(1,24);
    mem.p_mess->create_potR();
    mem.p_mess->zeroR(-frac.get_density_0());
    int howbig=2*mem.p_mess->total_memory;
    std::copy(mem.p_mess->potRS,mem.p_mess->potRS+howbig,mem.p_mess->potR);
    mem.p_mess->free_potRS();
    mem.p_mess->create_potC();
    mem.p_mess->fftw_real_to_complex();
    double g_c=group.get_force_const()/length_3;
    FileFFT << "going to power spectrum from periodic " << length_1 << "\n";
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
    for(int nx=mem.p_mess->start_x;nx < mem.p_mess->start_x+mem.p_mess->length_x;++nx)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nz=0;nz < length_c;++nz)
	      {
		int twit=mem.fftw_where(nx,ny,nz,length_1,length_c);
		double of_the_year=-g_c/(green_1[nz]+green_1[ny]+green_1[nx]);
		mem.p_mess->potC[twit][0]*=of_the_year;
		mem.p_mess->potC[twit][1]*=of_the_year;
	      }
	  }
      }
    mem.p_mess->fftw_complex_to_real();
    mem.p_mess->free_potC();
    mem.p_mess->create_potRS();
    howbig=2*mem.p_mess->total_memory;
    std::copy(mem.p_mess->potR,mem.p_mess->potR+howbig,mem.p_mess->potRS);
    mem.p_mess->free_potR();
    Full_Stop(mem,35);
    frac.timing(-1,24);
    slices_to_potf(mem,frac,0);
    frac.timing(1,24);
    FileFFT << "exiting periodic " << "\n";
    frac.timing(1,4);
  }
}

