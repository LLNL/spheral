#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group, Fractal& fractal,Misc& misc)
  {
    if(!fractal.amIaFftNode()) return;
    const int length_1=fractal.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    const double g_c=pow((double)length_1,-5)/8.0;
      //    const double g_c=1.0/(double)Misc::pow(length_1,2);
    size_t sizeR=sizeof(double);
    static vector <double> green(length_11*length_11*length_11);
    if(Fractal::first_time_solver)
      {
	double* greenR;
	fftw_complex* greenC;
	greenR=(double*) fftw_malloc(sizeR*2*length_11*length_2*length_2);
	greenC=(fftw_complex*) greenR;
	fftw_plan plan_green_cr=fftw_plan_dft_c2r_3d(length_1,length_1,length_1,greenC,greenR,FFTW_ESTIMATE);
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
		    greenR[fftw_where(n_x,n_y,n_z,length_2,length_2)]=dead_parrot;
		    greenR[fftw_where(n_x_inv,n_y,n_z,length_2,length_2)]=dead_parrot;
		    greenR[fftw_where(n_x,n_y_inv,n_z,length_2,length_2)]=dead_parrot;
		    greenR[fftw_where(n_x_inv,n_y_inv,n_z,length_2,length_2)]=dead_parrot;

		    greenR[fftw_where(n_x,n_y,n_z_inv,length_2,length_2)]=dead_parrot;
		    greenR[fftw_where(n_x_inv,n_y,n_z_inv,length_2,length_2)]=dead_parrot;
		    greenR[fftw_where(n_x,n_y_inv,n_z_inv,length_2,length_2)]=dead_parrot;
		    greenR[fftw_where(n_x_inv,n_y_inv,n_z_inv,length_2,length_2)]=dead_parrot;
		  }	      
	      }
	  }
	fftw_execute(plan_green_cr);
	for( int p_x=0;p_x<length_11;++p_x)
	  {
	    for( int p_y=0;p_y<length_11;++p_y)
	      {
		for( int p_z=0;p_z<length_11;++p_z)
		  {
		    //		    green(p_x,p_y,p_z)=real(green_c(p_x,p_y,p_z));
		    green[fftw_where(p_x,p_y,p_z,length_11,length_11)]=greenR[2*fftw_where(p_x,p_y,p_z,length_2,length_11)];
		  }
	      }
	  }
	fftw_free(greenR);
	greenR=0;
	greenC=0;
	fftw_destroy_plan(plan_green_cr);
	cout << "isolated what " << endl;
	Fractal::first_time_solver=false;
	return;
      }
    //
    cout << "isol 0 " << endl;
    double* potR;
    fftw_complex* potC;
    potR=(double*) fftw_malloc(sizeR*2*length_11*length_2*length_2);
    potC=(fftw_complex*) potR;
    static fftw_plan plan_rc;
    static fftw_plan plan_cr;
    plan_rc=fftw_plan_dft_r2c_3d(length_1,length_1,length_1,potR,potC,FFTW_ESTIMATE);
    plan_cr=fftw_plan_dft_c2r_3d(length_1,length_1,length_1,potC,potR,FFTW_ESTIMATE);
    for( int p_x=0;p_x<length_2;++p_x)
      {
	for( int p_y=0;p_y<length_2;++p_y)
	  {
	    for( int p_z=0;p_z<length_2;++p_z)
	      {
		potR[fftw_where(p_x,p_y,p_z,length_2,length_2)]=0.0;
	      }
	  }
      }
    for( int p_x=0;p_x<length_11;++p_x)
      {
	for( int p_y=0;p_y<length_11;++p_y)
	  {
	    for( int p_z=0;p_z<length_11;++p_z)
	      {
		int p_0=fftw_where(p_x,p_y,p_z,length_11,length_11);
		potR[fftw_where(p_x,p_y,p_z,length_2,length_2)]=group.list_points[p_0]->get_density_point();
	      }
	  }
      }
    cout << "isol a " << endl;
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
		int wherexyz=2*fftw_where(n_x,n_y,n_z,length_2,length_11);
		int whereabc=fftw_where(n_a,n_b,n_c,length_11,length_11);
		potR[wherexyz]*=green[whereabc];
		potR[wherexyz+1]*=green[whereabc];
		//		potC[fftw_where(n_x,n_y,n_z,length_2,length_11)]*=green[fftw_where(n_a,n_b,n_c,length_11,length_11)];
	      }
	  }
      }
    cout << "isol b " << endl;
    fftw_execute(plan_cr);
    cout << "isol c " << endl;
    for (int n_x=0;n_x <=length_1;++n_x)
      {
	for (int n_y=0;n_y <=length_1;++n_y)
	  {
	    for (int n_z=0;n_z <=length_1;++n_z)
	      {
		int n_p=fftw_where(n_x,n_y,n_z,length_11,length_11);
		group.list_points[n_p]->set_potential_point(potR[fftw_where(n_x,n_y,n_z,length_2,length_2)]);
	      }
	  }
      }
    fftw_free(potR);
    potR=0;
    potC=0;
  }
}
