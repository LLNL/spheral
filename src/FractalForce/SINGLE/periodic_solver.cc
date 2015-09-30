#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void periodic_solver(Group& group, Fractal_Memory& fractal_memory,Fractal& fractal)
  {
    cout << "entering periodic " << endl;
    const int length_1=fractal.get_grid_length();  
    const int length_c=(length_1+2)/2;
    const int length_2=length_1*length_1;
    const double pi=4.0*atan(1.0);
    group.set_grav_const(4.0*pi/(double)length_2);
    //
    size_t align=sizeof(Complex);
    vector <double> green_1(length_1+1);
    //
    for (int k=0;k <length_1+1;++k)
      {
	double aa=pi*(double)(k)/(double)(length_1);
	green_1[k]=1.0e-30+pow(2.0*sin(aa),2);
      }
    array3 <Complex> pot(length_1,length_1,length_c,align);
    rcfft3d Forward3  (pot);
    crfft3d Backward3 (pot);
    for(int k=0;k <length_1/2;++k)
      {
	for(int j=0;j <length_1;++j)
	  {
	    for(int i=0;i <length_1;++i)
	      {
		int ir=i+(j+2*k*length_1)*length_1;
		Point* p_point_a=group.list_points[ir];
		Point* p_point_b=group.list_points[ir+length_2];
		pot(i,j,k)=Complex(p_point_a->get_density_point(),p_point_b->get_density_point());
	      }
	  }
      }
    cout << "going to four_2 a " << endl;
    Forward3.fft(pot);
    cout << "going to power spectrum from periodic " << length_1 << endl;
    vector <double> variance_rho(length_1,0.0);
    vector <double> variance_pot(length_1,0.0);
    vector <double> variance_force(length_1,0.0);
    vector <double> variance_force_s(length_1,0.0);
    pot(0,0,0)=0.0;
    fractal.timing(-1,6);
    power_spectrum(pot,length_1,variance_rho,variance_pot,variance_force,variance_force_s,0,fractal.get_density_0(),fractal_memory.do_var);
    fractal.timing(1,6);
    double g_c=group.get_grav_const();
    cout << "going to four_2 b " << endl;
    for(int nz=0;nz < length_c;++nz)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nx=0;nx < length_1;++nx)
	      {
		pot(nx,ny,nz)*=-g_c/(green_1[nz]+green_1[ny]+green_1[nx]);
	      }
	  }
      }
    cout << "going to four_2 c " << endl;
    Backward3.fftNormalized(pot);
    for(int nz=0;nz < length_1/2;++nz)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nx=0;nx < length_1;++nx)
	      {
		int ir=nx+(ny+2*nz*length_1)*length_1;
		Point* p_point_a=group.list_points[ir];
		p_point_a->set_potential_point(real(pot(nx,ny,nz)));
		Point* p_point_b=group.list_points[ir+length_2];
		p_point_b->set_potential_point(imag(pot(nx,ny,nz)));
	      }
	  }
      }
    cout << "exiting periodic " << endl;
  }
}
