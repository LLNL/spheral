
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void max_predict(Fractal_Memory& mem,Fractal& fractal,vector <double>& shear_force,double& r)
  {
    vector <double>shear(6);
    vector <double> displ(6);
    double vol_crash=1.0/mem.density_crash;
    double vol_old=1.0;
    for(int i=0;i<6;i++)
      shear[i]=shear_force[i]*fractal.omega_fraction;
    double min_vol=1.0;
    r=min_vol;
    double dr=fractal.rad[1]-fractal.rad[0];
    for(int i=0;i<=100;i++)
      {
	double haha=fractal.grow[i];
	displ[0]=1.0+haha*shear[0];
	displ[1]=haha*shear[1];
	displ[2]=haha*shear[2];
	displ[3]=1.0+haha*shear[3];
	displ[4]=haha*shear[4];
	displ[5]=1.0+haha*shear[5];
	double vol=det(displ);
	if(vol < min_vol)
	  {
	    min_vol=vol;
	    r=fractal.grow[i];
	    if(min_vol < vol_crash)
	      {
		if(i > 0)
		  dr=fractal.rad[i]-fractal.rad[i-1];
		double dvol=(vol-vol_old)/dr;
		double dx=-(vol-vol_crash)/dvol;
		//		cerr << "maxx " << vol << " "  << vol_old << " "  << vol_crash << " "  << dr << " "  << dx << "\n";
		r=-(fractal.rad[i]+dx);
		assert(dx/dr > -1.0 && dx <= 0.0);
		//		cerr << "maxx " << r << " " << dx << "\n";
		return;
	      }
	  }
	vol_old=vol;
      }
  }
  template <class T> T det(const vector <T>& m)
  {
    return m[0]*m[3]*m[5]+2.0*m[1]*m[4]*m[2]-
      m[2]*m[2]*m[3]-m[0]*m[4]*m[4]-m[1]*m[1]*m[5];
  }
}
