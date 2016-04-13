#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void periodic_solver(Group& group, Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFractal=frac.p_file->DUMPS;
    //    ofstream& FileFractal=frac.p_file->FileFractal;
    FileFractal << "entering periodic " << "\n";
    const int length_1=frac.get_grid_length();  
    const int length_c=(length_1+2)/2;
    const int length_2=length_1*length_1;
    const double length_3=static_cast<double>(length_1*length_2);
    const double pi=4.0*atan(1.0);
    group.set_force_const(4.0*pi/(double)length_2);
    //
    size_t sizeR=sizeof(double);
    size_t sizeC=sizeof(fftw_complex);
    double* potR;
    fftw_complex* potC;
    potR=(double*) fftw_malloc(sizeR*2*length_c*length_2);
    potC=(fftw_complex*) fftw_malloc(sizeC*length_c*length_2);
    fftw_plan plan_rc=fftw_plan_dft_r2c_3d(length_1,length_1,length_1,potR,potC,FFTW_ESTIMATE);
    fftw_plan plan_cr=fftw_plan_dft_c2r_3d(length_1,length_1,length_1,potC,potR,FFTW_ESTIMATE);
    vector <double> green_1(length_1+1);
    //
    for (int k=0;k <length_1+1;++k)
      {
	double aa=pi*(double)(k)/(double)(length_1);
	green_1[k]=1.0e-30+pow(2.0*sin(aa),2);
      }
    double summ0=0.0;
    double summ1=0.0;
    double summ2=0.0;
    for(int i=0;i <length_1;++i)
      {
	for(int j=0;j <length_1;++j)
	  {
	    for(int k=0;k <length_1;++k)
	      {
		int ir=frac.where_1(i,j,k);
		Point* p_point_a=group.list_points[ir];
		potR[mem.fftw_where(i,j,k,length_1,length_1)]=p_point_a->get_density_point();
		summ0++;
		summ1+=potR[mem.fftw_where(i,j,k,length_1,length_1)];
		summ2+=pow(potR[mem.fftw_where(i,j,k,length_1,length_1)],2);
		/*
		FileFractal << "dens " ;
		FileFractal << i << " " ;
		FileFractal << j << " " ;
		FileFractal << k << " " ;
		FileFractal << mem.fftw_where(i,j,k,length_1,length_1) << " " ;
		FileFractal << potR[mem.fftw_where(i,j,k,length_1,length_1)] << "\n";
		*/
	      }
	  }
      }
    FileFractal << "going to four_2 a " << "\n";
    //    assert(0);
    summ1/=summ0;
    summ2/=summ0;
    summ2=sqrt(summ2);
    FileFractal << "summ " << summ0 << " " << summ1 << " " << summ2 << "\n";
    fftw_execute(plan_rc);
    double g_c=group.get_force_const()/length_3;
    FileFractal << "going to power spectrum from periodic " << length_1 << "\n";
    vector <double> variance_rho(length_1,0.0);
    vector <double> variance_pot(length_1,0.0);
    vector <double> variance_force(length_1,0.0);
    vector <double> variance_force_s(length_1,0.0);
    potC[0][0]=0.0;
    frac.timing(-1,6);
    power_spectrum(potC,length_1,variance_rho,variance_pot,variance_force,variance_force_s,0,frac.get_density_0(),mem.do_var,mem);
    frac.timing(1,6);
    FileFractal << "going to four_2 b " << "\n";
    for(int nx=0;nx < length_1;++nx)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nz=0;nz < length_c;++nz)
	      {
		int twit=mem.fftw_where(nx,ny,nz,length_1,length_c);
		double of_the_year=-g_c/(green_1[nz]+green_1[ny]+green_1[nx]);
		potC[twit][0]*=of_the_year;
		potC[twit][1]*=of_the_year;
	      }
	  }
      }
    FileFractal << "going to four_2 c " << "\n";
    fftw_execute(plan_cr);
    FileFractal << "level zero potential a " << "\n";
    for(int nx=0;nx < length_1;++nx)
      {
	for(int ny=0;ny < length_1;++ny)
	  {
	    for(int nz=0;nz < length_1;++nz)
	      {
		int ir=frac.where_1(nx,ny,nz);
		Point* p_point_a=group.list_points[ir];
		p_point_a->set_potential_point(potR[mem.fftw_where(nx,ny,nz,length_1,length_1)]);
	      }
	  }
      }
    FileFractal << "level zero potential b " << "\n";
    fftw_free(potR);
    fftw_free(potC);
    potR=0;
    potC=0;
    fftw_destroy_plan(plan_rc);
    fftw_destroy_plan(plan_cr);
    Fractal::first_time_solver=false;
    FileFractal << "exiting periodic " << "\n";
  }
}
