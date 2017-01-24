#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFractal=frac.p_file->FileFractal;
    const int length_1=frac.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    const double g_c=pow((double)length_1,-5)/8.0;
    //    const double g_c=1.0/(double)Misc::pow(length_1,2);
    size_t sizeR=sizeof(double);
    size_t sizeC=sizeof(fftw_complex);
    static vector <double> green(length_11*length_11*length_11);
    if(Fractal::first_time_solver)
      {
	double* greenR;
	fftw_complex* greenC;
	greenR=(double*) fftw_malloc(sizeR*2*length_11*length_2*length_2);
	greenC=(fftw_complex*) fftw_malloc(sizeC*length_11*length_2*length_2);
	fftw_plan plan_green_rc=fftw_plan_dft_r2c_3d(length_2,length_2,length_2,greenR,greenC,FFTW_ESTIMATE);
	for(int n_x=0;n_x < length_11;++n_x)
	  {
	    int n_x_inv=n_x;
	    if(n_x > 0) n_x_inv=length_2-n_x;
	    double x_2=(double)(n_x*n_x);
	    for(int n_y=0;n_y < length_11;++n_y)
	      {
		int n_y_inv=n_y;
		if(n_y > 0) n_y_inv=length_2-n_y;
		double y_2=(double)(n_y*n_y);
		for( int n_z=0;n_z<length_11;++n_z)
		  {
		    int n_z_inv=n_z;
		    if(n_z > 0) n_z_inv=length_2-n_z;
		    double z_2=(double)(n_z*n_z)+0.25;
		    double r_2=z_2+y_2+x_2;
		    double dead_parrot=-g_c/sqrt(r_2);
		    greenR[mem.fftw_where(n_x,n_y,n_z,length_2,length_2)]=dead_parrot;
		    greenR[mem.fftw_where(n_x_inv,n_y,n_z,length_2,length_2)]=dead_parrot;
		    greenR[mem.fftw_where(n_x,n_y_inv,n_z,length_2,length_2)]=dead_parrot;
		    greenR[mem.fftw_where(n_x_inv,n_y_inv,n_z,length_2,length_2)]=dead_parrot;

		    greenR[mem.fftw_where(n_x,n_y,n_z_inv,length_2,length_2)]=dead_parrot;
		    greenR[mem.fftw_where(n_x_inv,n_y,n_z_inv,length_2,length_2)]=dead_parrot;
		    greenR[mem.fftw_where(n_x,n_y_inv,n_z_inv,length_2,length_2)]=dead_parrot;
		    greenR[mem.fftw_where(n_x_inv,n_y_inv,n_z_inv,length_2,length_2)]=dead_parrot;
		  }	      
	      }
	  }
	fftw_execute(plan_green_rc);
	for( int p_x=0;p_x<length_11;++p_x)
	  {
	    for( int p_y=0;p_y<length_11;++p_y)
	      {
		for( int p_z=0;p_z<length_11;++p_z)
		  {
		    green[mem.fftw_where(p_x,p_y,p_z,length_11,length_11)]=greenC[mem.fftw_where(p_x,p_y,p_z,length_2,length_11)][0];
		  }
	      }
	  }
	fftw_free(greenR);
	fftw_free(greenC);
	greenR=0;
	greenC=0;
	fftw_destroy_plan(plan_green_rc);
	FileFractal << "isolated what " << "\n";
	Fractal::first_time_solver=false;
	return;
      }
    //
    FileFractal << "isol 0 " << "\n";
    double* potR;
    fftw_complex* potC;
    potR=(double*) fftw_malloc(sizeR*2*length_11*length_2*length_2);
    potC=(fftw_complex*) fftw_malloc(sizeR*2*length_11*length_2*length_2);
    fftw_plan plan_rc=fftw_plan_dft_r2c_3d(length_2,length_2,length_2,potR,potC,FFTW_ESTIMATE);
    fftw_plan plan_cr=fftw_plan_dft_c2r_3d(length_2,length_2,length_2,potC,potR,FFTW_ESTIMATE);

    for( int p_x=0;p_x<length_2;++p_x)
      {
	for( int p_y=0;p_y<length_2;++p_y)
	  {
	    for( int p_z=0;p_z<length_2;++p_z)
	      {
		potR[mem.fftw_where(p_x,p_y,p_z,length_2,length_2)]=0.0;
	      }
	  }
      }
    for( int p_x=0;p_x<length_11;++p_x)
      {
	for( int p_y=0;p_y<length_11;++p_y)
	  {
	    for( int p_z=0;p_z<length_11;++p_z)
	      {
		//		int p_0=mem.fftw_where(p_x,p_y,p_z,length_11,length_11);
		int p_0=frac.where_1(p_x,p_y,p_z);
		potR[mem.fftw_where(p_x,p_y,p_z,length_2,length_2)]=group.list_points[p_0]->get_density_point();
	      }
	  }
      }
    FileFractal << "isol a " << "\n";
    fftw_execute(plan_rc);
    //
    for(int n_x=0;n_x < length_2;++n_x)
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
		potC[wherexyz][0]*=green[whereabc];
		potC[wherexyz][1]*=green[whereabc];
	      }
	  }
      }
    FileFractal << "isol b " << "\n";
    fftw_execute(plan_cr);
    FileFractal << "isol c " << "\n";
    for (int n_x=0;n_x < length_11;++n_x)
      {
	for (int n_y=0;n_y < length_11;++n_y)
	  {
	    for (int n_z=0;n_z < length_11;++n_z)
	      {
		//		int n_p=mem.fftw_where(n_x,n_y,n_z,length_11,length_11);
		int n_p=frac.where_1(n_x,n_y,n_z);
		group.list_points[n_p]->set_potential_point(potR[mem.fftw_where(n_x,n_y,n_z,length_2,length_2)]);
	      }
	  }
      }
    fftw_free(potR);
    fftw_free(potC);
    potR=0;
    potC=0;
    fftw_destroy_plan(plan_rc);
    fftw_destroy_plan(plan_cr);
  }
}
