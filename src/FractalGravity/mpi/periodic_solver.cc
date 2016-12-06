#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void periodic_solver(Group& group, Fractal_Memory& fractal_memory,Fractal& fractal)
  {
    if(!fractal.amIaFftNode()) return;
    cout << "entering periodic " << endl;
    const int length_1=fractal.get_grid_length();  
    const int length_c=(length_1+2)/2;
    const int length_2=length_1*length_1;
    const double length_3=(double)(length_1*length_2);
    const double pi=4.0*atan(1.0);
    group.set_grav_const(4.0*pi/(double)length_2);
    //
    size_t sizeR=sizeof(double);
    double* potR;
    fftw_complex* potC;
    potR=(double*) fftw_malloc(sizeR*2*length_c*length_2);
    potC=(fftw_complex*) potR;
    static fftw_plan plan_rc;
    static fftw_plan plan_cr;
    if(Fractal::first_time_solver)
      {
	plan_rc=fftw_plan_dft_r2c_3d(length_1,length_1,length_1,potR,potC,FFTW_ESTIMATE);
	plan_cr=fftw_plan_dft_c2r_3d(length_1,length_1,length_1,potC,potR,FFTW_ESTIMATE);
      }
    //
    vector <double> green_1(length_1+1);
    vector <int>Box(6);
    fractal.getBox(Box);
    vector <int>BoxReal(6);
    fractal.getBoxReal(BoxReal);
    //
    for (int k=0;k <length_1+1;++k)
      {
	double aa=pi*(double)(k)/(double)(length_1);
	green_1[k]=1.0e-30+pow(2.0*sin(aa),2);
      }
    for(int i=0;i <length_1;++i)
      {
	for(int j=0;j <length_1;++j)
	  {
	    for(int k=0;k <length_1;++k)
	      {
		int ir=where_1(i,j,k,Box);
		Point* p_point_a=group.list_points[ir];
		potR[fftw_where(i,j,k,length_1,length_1)]=p_point_a->get_density_point();
	      }
	  }
      }
    cout << "going to four_2 a " << endl;
    fftw_execute(plan_rc);
    double g_c=group.get_grav_const()/length_3;
    cout << "going to power spectrum from periodic " << length_1 << endl;
    vector <double> variance_rho(length_1,0.0);
    vector <double> variance_pot(length_1,0.0);
    vector <double> variance_force(length_1,0.0);
    vector <double> variance_force_s(length_1,0.0);
    potR[0]=0.0;
    potR[1]=0.0;
    fractal.timing(-1,6);
    power_spectrum(potR,length_1,variance_rho,variance_pot,variance_force,variance_force_s,0,fractal.get_density_0(),fractal_memory.do_var);
    fractal.timing(1,6);
    cout << "going to four_2 b " << endl;
    for(int nx=0;nx < length_1;++nx)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nz=0;nz < length_c;++nz)
	      {
		int twit=2*fftw_where(nx,ny,nz,length_1,length_c);
		double of_the_year=-g_c/(green_1[nz]+green_1[ny]+green_1[nx]);
		potR[twit]*=of_the_year;
		potR[twit+1]*=of_the_year;
		//		potC[fftw_where(nx,ny,nz,length_1,length_c)]*=-g_c/(green_1[nz]+green_1[ny]+green_1[nx]);
		//		pot[nx][ny][nz]*=-g_c/(green_1[nz]+green_1[ny]+green_1[nx]);
	      }
	  }
      }
    cout << "going to four_2 c " << endl;
    fftw_execute(plan_cr);
    cout << "level zero potential a " << endl;
    for(int nx=0;nx < length_1;++nx)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nz=0;nz < length_1;++nz)
	      {
		int ir=where_1(nx,ny,nz,Box);
		Point* p_point_a=group.list_points[ir];
		p_point_a->set_potential_point(potR[fftw_where(nx,ny,nz,length_1,length_1)]);
	      }
	  }
      }
    cout << "level zero potential b " << endl;
    for(int nx=Box[0];nx <= Box[1];++nx)
      {
	for(int ny=Box[2];ny <= Box[3];++ny)
	  {
	    for(int nz=Box[4];nz <= Box[5];++nz)
	      {
		int ir=where_1(nx,ny,nz,Box);
		Point* p_point=group.list_points[ir];
		if(p_point->get_buffer_point())
		  {
		    int nx1=(nx+length_1) % length_1;
		    int ny1=(ny+length_1) % length_1;
		    int nz1=(nz+length_1) % length_1;
		    int ir1=where_1(nx1,ny1,nz1,Box);
		    Point* p_point_1=group.list_points[ir1];
		    p_point->set_potential_point(p_point_1->get_potential_point());
		    p_point_1->dumpp();
		    p_point->dumpp();
		  }
	      }
	  }
      }
    cout << "level zero potential c " << endl;
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	point.dumpp();
      }
    fftw_free(potR);
    potR=0;
    potC=0;
    Fractal::first_time_solver=false;
    cout << "exiting periodic " << endl;
  }
  int fftw_where(const int& i,const int& j,const int& k,const int& la,const int& lb)
  {
    return k+(j+i*la)*lb;
  }
}
