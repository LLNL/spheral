#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    frac.timing(-1,5);
    static bool printit=true;
    Fractal::first_time_solver=false;
    printit=false;
    ofstream& FileFFT=frac.p_file->DUMPS;
    FileFFT << "entering isolated " << "\n";
    int length_1=frac.get_grid_length();
    int length_11=length_1+1;
    int length_2=2*length_1;
    int length_22=length_11*2;
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
		if(printit && abs(mem.p_mess->potR[wherexyz]) > 1.0e-5) FileFFT << " SP0 " << wherexyz << " " << nxyz << " " << nx << " " << ny << " " << nz << " " << mem.p_mess->potR[wherexyz] << "\n";
		nxyz++;
	      }
	  }
      }
    //
    //    double close_your_eyes=1.0/static_cast<double>(length_2*length_2*length_2);
    //
    mem.p_mess->free_potRS();
    mem.p_mess->create_potC();
    for(int nx=mem.p_mess->start_x;nx < mem.p_mess->start_x+mem.p_mess->length_x;++nx)
      {
	for(int ny=0;ny < length_2;++ny)
	  {
	    for(int nz=0;nz < length_22;++nz)
	      {
		int wherexyz=mem.fftw_where(nx,ny,nz,length_2,length_22);
		if(printit && abs(mem.p_mess->potR[wherexyz]) > 1.0e-5) FileFFT << " SPM " << wherexyz << " " << nx << " " << ny << " " << nz << " " << mem.p_mess->potR[wherexyz] << "\n";
	      }
	  }
      }
    mem.p_mess->fftw_real_to_complex();
    if(mem.p_mess->FractalRank == 0)
      FileFFT << " zero power " << mem.p_mess->potC[0][0] << " " << mem.p_mess->potC[0][1] << "\n";
    //    double cheese_shop=1.0; /////////////
    for(int nx=mem.p_mess->start_x;nx < mem.p_mess->start_x+mem.p_mess->length_x;++nx)
      {
	for(int ny=0;ny < length_2;++ny)
	  {
	    int nb=min(ny,length_2-ny);
	    for(int nz=0;nz < length_11;++nz)
	      {
		int nc=min(nz,length_2-nz);
		int wherexyz=mem.fftw_where(nx,ny,nz,length_2,length_11);
		int whereabc=mem.fftw_where(nx,nb,nc,length_11,length_11);
		assert(wherexyz >= 0);
		assert(whereabc >= 0);
		mem.p_mess->potC[wherexyz][0]*=mem.p_mess->green[whereabc];
		mem.p_mess->potC[wherexyz][1]*=mem.p_mess->green[whereabc];
		//		mem.p_mess->potC[wherexyz][0]*=cheese_shop;   ///////////////
		//		  mem.p_mess->potC[wherexyz][1]*=cheese_shop;   /////////////
		double sum=sqrt(pow(mem.p_mess->potC[wherexyz][0],2)+pow(mem.p_mess->potC[wherexyz][1],2));
		if(printit) FileFFT << " SPC " << wherexyz << " " << nx << " " << ny << " " << nz << " " << mem.p_mess->potC[wherexyz][0] << " " << mem.p_mess->potC[wherexyz][1] << " " << sum <<"\n";
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
		if(printit && abs(mem.p_mess->potR[nxyz]) > 1.0e-5) FileFFT << " SP1 " << wherexyz << " " << nxyz << " " << nx << " " << ny << " " << nz << " " << mem.p_mess->potRS[nxyz] << " " << mem.p_mess->potR[wherexyz] << "\n";
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
    printit=false;
    frac.timing(1,5);
  }
}
 
