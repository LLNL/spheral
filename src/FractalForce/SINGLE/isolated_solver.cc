#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void isolated_solver(Group& group, Fractal& fractal,Misc& misc)
  {
    const int length_ratio=fractal.get_length_ratio();
    const int length_1=fractal.get_grid_length();
    const int length_11=length_1+1;
    const int length_2=2*length_1;
    int length_z=length_1;
    int length_2z=2*length_z;
    if(length_ratio != 1)
      {
	length_z/=length_ratio;
	length_2z/=length_ratio;      
      }
    const int length_1z=length_z+1;
    const double g_c=1.0/(double)Misc::pow(length_1,2);
    size_t align=sizeof(Complex);
    static array3 <double> green(length_11,length_11,length_1z);
    if(Fractal::first_time_solver)
      {
	array3 <double> rho(length_2,length_2,length_2z,align) ;
	array3 <Complex> green_c(length_2,length_2,(length_2z+2)/2,align) ;
	for( int n_z=0;n_z<length_1z;++n_z)
	  {
	    int n_z_inv=n_z;
	    if(n_z > 0) n_z_inv=length_2z-n_z;
	    double z_2=(double)(n_z*n_z)+0.25;
	    for(int n_y=0;n_y < length_11;++n_y)
	      {
		int n_y_inv=n_y;
		if(n_y > 0) n_y_inv=length_2-n_y;
		double y_2=(double)(n_y*n_y);
		for(int n_x=0;n_x < length_11;++n_x)
		  {
		    int n_x_inv=n_x;
		    if(n_x > 0) n_x_inv=length_2-n_x;
		    double r_2=z_2+y_2+(double)(n_x*n_x);
		    rho(n_x,n_y,n_z)=-g_c/sqrt(r_2);
		    rho(n_x_inv,n_y,n_z)=rho(n_x,n_y,n_z);
		    rho(n_x,n_y_inv,n_z)=rho(n_x,n_y,n_z);
		    rho(n_x_inv,n_y_inv,n_z)=rho(n_x,n_y,n_z);
		    rho(n_x,n_y,n_z_inv)=rho(n_x,n_y,n_z);
		    rho(n_x_inv,n_y,n_z_inv)=rho(n_x,n_y,n_z);
		    rho(n_x,n_y_inv,n_z_inv)=rho(n_x,n_y,n_z);
		    rho(n_x_inv,n_y_inv,n_z_inv)=rho(n_x,n_y,n_z);
		  }	      
	      }
	  }
	rcfft3d Forward3 (rho,green_c);
	Forward3.fft(rho,green_c);
	for( int p_z=0;p_z<length_1z;++p_z)
	  {
	    for( int p_y=0;p_y<length_11;++p_y)
	      {
		for( int p_x=0;p_x<length_11;++p_x)
		  {
		    green(p_x,p_y,p_z)=real(green_c(p_x,p_y,p_z));
		  }
	      }
	  }
	cout << "isolated what " << endl;
	Fractal::first_time_solver=false;
	return;
      }
    //
    cout << "isol 0 " << endl;
    array3 <double> rho(length_2,length_2,length_2z,align);
    array3 <Complex> pot(length_2,length_2,length_1z,align);
    for( int p_z=0;p_z<length_2z;++p_z)
      {
	for( int p_y=0;p_y<length_2;++p_y)
	  {
	    for( int p_x=0;p_x<length_2;++p_x)
	      {
		rho(p_x,p_y,p_z)=0.0;
	      }
	  }
      }
    for( int p_z=0;p_z<length_1z;++p_z)
      {
	for( int p_y=0;p_y<length_11;++p_y)
	  {
	    for( int p_x=0;p_x<length_11;++p_x)
	      {
		int p_0=Misc::nr(p_x,p_y,p_z,length_11);
		rho(p_x,p_y,p_z)=group.list_points[p_0]->get_density_point();
	      }
	  }
      }
    cout << "isol a " << endl;
    rcfft3d Forward3 (rho,pot);
    Forward3.fft(rho,pot);
    //
    for(int n_z=0;n_z < length_1z;++n_z)
      {
	int n_c=min(n_z,length_2z-n_z);
	for(int n_y=0;n_y < length_2;++n_y)
	  {
	    int n_b=min(n_y,length_2-n_y);
	    for(int n_x=0;n_x < length_2;++n_x)
	      {
		int n_a=min(n_x,length_2-n_x);
		pot(n_x,n_y,n_z)=pot(n_x,n_y,n_z)*green(n_a,n_b,n_c);
	      }
	  }
      }
    cout << "isol b " << endl;
    crfft3d Backward3 (pot,rho);
    Backward3.fftNormalized(pot,rho);
    cout << "isol c " << endl;
    for (int n_z=0;n_z <=length_z;++n_z)
      {
	for (int n_y=0;n_y <=length_1;++n_y)
	  {
	    for (int n_x=0;n_x <=length_1;++n_x)
	      {
		int n_p=Misc::nr(n_x,n_y,n_z,length_11);
		group.list_points[n_p]->set_potential_point(rho(n_x,n_y,n_z));
	      }
	  }
      }
  }
}
