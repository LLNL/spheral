#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <class M, class F>  void make_particles(M& mem,F& frac,int& count,const double& m,const bool&)
  {
    ofstream& FileFractal=mem.p_file->DUMPS;
    //    ofstream& FileFractal=mem.p_file->FileFractal;
    //    ofstream& FilePos=mem.p_file->FilePos;
    FileFractal << " particle a m= " << m << "\n";
    double rand_max=(double)RAND_MAX;
    //    double rmax=1.0e-5;
    double rmax=0.3;
    double x_off=0.405;
    double y_off=0.46;
    double z_off=0.4247;
    double slope=-0.5;
    //
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    //
    int nparts=frac.get_number_particles();

    //    nparts=2;

    double rmax2=rmax*rmax;
    double total_mass=m*static_cast<double>(nparts);
    double sigma2=total_mass/(2.0*rmax);
    //    sigma2=0.0; //*************//
    try
      {
	Particle* particles=new Particle[nparts];
      }
    catch(bad_alloc& ba)
      {
	cerr << " bad galaxy allocation " << nparts << " " << ba.what() << endl;
	exit(0);
      }
    frac.particle_list.resize(nparts);
    FileFractal << " parta " << nparts << "\n";
    frac.set_number_particles(nparts);
    FileFractal << " partb " << nparts << "\n";
    vector <double>vel(3,0.0);
    for(int n=0;n<nparts;++n)
      {
	frac.particle_list[n]=&particles[n];
	particles[n].space_resize(6);
	particles[n].field_resize(4);
	double r1=Fractal::my_rand(rand_max);
	r1=pow(r1,expo);
	double ctheta=2.0*Fractal::my_rand(rand_max)-1.0;
	double phi=twopi*Fractal::my_rand(rand_max);
	double r=rmax*r1;
	double z=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double x=r*stheta*cos(phi)+x_off;
	double y=r*stheta*sin(phi)+y_off;
	frac.particle_list[n]->set_pos(x,y,z);
	double sig2=sigma2*(1.0-r*r/rmax2);
	double sig=sqrt(sig2);
	double vel2=sig*sqrt(-2.0*log(Fractal::my_rand(rand_max)));
	phi=twopi*Fractal::my_rand(rand_max);
	vel[0]=vel2*cos(phi);
	vel[1]=vel2*sin(phi);
	vel2=sig*sqrt(-2.0*log(Fractal::my_rand(rand_max)));
	phi=twopi*Fractal::my_rand(rand_max);
	vel[2]=vel2*cos(phi);
	/*
	if(mem.p_mess->FractalRank == 1)
	  vel[0]=0.6;
	else
	  vel[0]=-0.6;
	*/
	frac.particle_list[n]->set_vel(vel);
	frac.particle_list[n]->set_mass(m);
	//	FilePos << " init " << n << "\t" << r << "\t" << sig << "\t" << x << "\t" << y << "\t" << z << "\t";
	//	FilePos << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t" << "\n";
      }
    count=nparts;
    FileFractal << " particle b m= " << m << "\n";
  }
}
namespace FractalSpace
{
  template void make_particles(Fractal_Memory& mem,Fractal& frac,int& count,const double& m,const bool& crash);
}
