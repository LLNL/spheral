#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFFT=frac.p_file->DUMPS;
    FileFFT << "entering isolated " << "\n";
    const int length_1=frac.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    int length_22=length_11*2;
    if(Fractal::first_time_solver)
      {
	Fractal::first_time_solver=false;
	return;
      }
    Full_Stop(mem,34);
    frac.timing(-1,24);
    dens_to_slices(group,mem,frac);
    frac.timing(1,24);
    mem.p_mess->create_potR();
    mem.p_mess->zeroR();
    int nxyz=0;
    for(int nx=mem.p_mess->start_x;nx < mem.p_mess->start_x+mem.p_mess->length_x;nx++)
      {
	for(int ny=0;ny<length_11;ny++)
	  {
	    for(int nz=0;nz<length_11;nz++)
	      {
		int wherexyz=mem.fftw_where(nx,ny,nz,length_2,length_22);
		mem.p_mess->potR[wherexyz]=mem.p_mess->potRS[nxyz];
		//		FileFFT << " SP0 " << wherexyz << " " << nxyz << " " << mem.p_mess->potRS[nxyz] << "\n";
		nxyz++;
	      }
	  }
      }
    //
    //    double close_your_eyes=1.0/static_cast<double>(length_2*length_2*length_2);
    //
    mem.p_mess->free_potRS();
    mem.p_mess->create_potC();
    mem.p_mess->fftw_real_to_complex();
    if(mem.p_mess->FractalRank == 0)
      FileFFT << " zero power " << mem.p_mess->potC[0][0] << " " << mem.p_mess->potC[0][1] << "\n";
    //
    for(int n_x=mem.p_mess->start_x;n_x < mem.p_mess->start_x+mem.p_mess->length_x;++n_x)
      {
	for(int n_y=0;n_y < length_2;++n_y)
	  {
	    int n_b=min(n_y,length_2-n_y);
	    for(int n_z=0;n_z < length_11;++n_z)
	      {
		int n_c=min(n_z,length_2-n_z);
		int wherexyz=mem.fftw_where(n_x,n_y,n_z,length_2,length_11);
		int whereabc=mem.fftw_where(n_x,n_b,n_c,length_11,length_11);
		assert(wherexyz >= 0);
		assert(whereabc >= 0);
		mem.p_mess->potC[wherexyz][0]*=mem.p_mess->green[whereabc];
		mem.p_mess->potC[wherexyz][1]*=mem.p_mess->green[whereabc];
		//		mem.p_mess->potC[wherexyz][0]*=close_your_eyes;
		//		mem.p_mess->potC[wherexyz][1]*=close_your_eyes;
	      }
	  }
      }
    mem.p_mess->fftw_complex_to_real();
    mem.p_mess->free_potC();
    mem.p_mess->create_potRS();
    nxyz=0;
    for(int nx=mem.p_mess->start_x;nx<mem.p_mess->start_x+mem.p_mess->length_x;nx++)
      {
	for(int ny=0;ny<length_11;ny++)
	  {
	    for(int nz=0;nz<length_11;nz++)
	      {
		int wherexyz=mem.fftw_where(nx,ny,nz,length_2,length_22);
		mem.p_mess->potRS[nxyz]=mem.p_mess->potR[wherexyz];
		nxyz++;
	      }
	  }
      }
    mem.p_mess->free_potR();
    Full_Stop(mem,35);
    frac.timing(-1,24);
    slices_to_potf(mem,frac,0);
    frac.timing(1,24);
    FileFFT << "exiting isolated " << "\n";
  }
}
