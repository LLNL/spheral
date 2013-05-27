#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFFT=frac.p_file->FileFFT;
    FileFFT << "entering isolated " << endl;
    const int length_1=frac.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    if(Fractal::first_time_solver)
      {
	Fractal::first_time_solver=false;
	return;
      }
    FileFFT << "isol 0 " << endl;
    mem.p_mess->create_potRC();
    mem.p_mess->zeroR();
    frac.timing(-1,24);
    dens_to_slices(group,mem,frac);
    frac.timing(1,24);
    FileFFT << "isol a " << endl;
    mem.p_mess->fftw_real_to_complex();
    if(mem.p_mess->FractalRank == 0)
      FileFFT << " zero power " << mem.p_mess->potC[0][0] << " " << mem.p_mess->potC[0][1] << endl;
    FileFFT << "isol aa " << endl;
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
	      }
	  }
      }
    FileFFT << "isol b " << endl;
    mem.p_mess->fftw_complex_to_real();
    FileFFT << "isol c " << endl;
    frac.timing(1,5);
    frac.timing(-1,24);
    slices_to_potf(group,mem,frac);
    frac.timing(1,24);
    frac.timing(-1,5);
    FileFFT << "isol d " << endl;
    mem.p_mess->free_potRC();
    FileFFT << "exiting isolated " << endl;
  }
}
