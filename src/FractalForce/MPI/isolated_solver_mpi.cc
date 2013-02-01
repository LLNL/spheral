#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFractal=frac.p_file->FileFractal;
    FileFractal << "entering isolated " << endl;
    const int length_1=frac.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    if(Fractal::first_time_solver)
      {
	Fractal::first_time_solver=false;
	return;
      }
    FileFractal << "isol 0 " << endl;
    mem.p_mess->create_potRC();
    mem.p_mess->zeroR();
    frac.timing(-1,24);
    dens_to_slices(group,mem,frac);
    frac.timing(1,24);
    FileFractal << "isol a " << endl;
    //    fftw_execute(mem.p_mess->plan_rc);
    //    mem.p_mess->dumpR(FileFractal,length_2,true);
    mem.p_mess->fftw_real_to_complex();
    //    fftw_mpi_execute_dft_r2c(mem.p_mess->plan_rc,mem.p_mess->potR,mem.p_mess->potC);
    if(mem.p_mess->FractalRank == 0)
      FileFractal << " zero power " << mem.p_mess->potC[0][0] << " " << mem.p_mess->potC[0][1] << endl;
    FileFractal << "isol aa " << endl;
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
		//		FileFractal << " where isol " << n_x << " "  << n_y << " "  << n_z << " " ;
		//		FileFractal << n_b << " "  << n_c << " " ;
		//		FileFractal << wherexyz << " "  << whereabc << " ";
		assert(wherexyz >= 0);
		assert(whereabc >= 0);
		//		mem.p_mess->potC[wherexyz][0]=262144.0*mem.p_mess->green[whereabc];
		//		mem.p_mess->potC[wherexyz][1]=0.0;
		mem.p_mess->potC[wherexyz][0]*=mem.p_mess->green[whereabc];
		mem.p_mess->potC[wherexyz][1]*=mem.p_mess->green[whereabc];
		//		FileFractal << mem.p_mess->green[whereabc] << " " << mem.p_mess->potC[wherexyz][0] << " " << mem.p_mess->potC[wherexyz][1] << endl;
	      }
	  }
      }
    FileFractal << "isol b " << endl;
    //fftw_execute(mem.p_mess->plan_cr);
    mem.p_mess->fftw_complex_to_real();
    //    fftw_mpi_execute_dft_c2r(mem.p_mess->plan_cr,mem.p_mess->potC,mem.p_mess->potR);
    FileFractal << "isol c " << endl;
    int length_22=length_2+2;
    /*
    for(int n_x=mem.p_mess->start_x;n_x < mem.p_mess->start_x+mem.p_mess->length_x;++n_x)
      {
	for(int n_y=0;n_y < length_11;++n_y)
	  {
	    for(int n_z=0;n_z < length_11;++n_z)
	      {
		int wherexyz=mem.fftw_where(n_x,n_y,n_z,length_2,length_22);
		mem.p_file->FileForce << n_x << " "  << n_y << " "  << n_z << " " << wherexyz << " " << mem.p_mess->potR[wherexyz] << endl;
	      }
	  }
      }
    */
    frac.timing(1,5);
    frac.timing(-1,24);
    slices_to_potf(group,mem,frac);
    frac.timing(1,24);
    frac.timing(-1,5);
    FileFractal << "isol d " << endl;
    mem.p_mess->free_potRC();
    FileFractal << "exiting isolated " << endl;
  }
}
