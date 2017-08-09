#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_me_a_galaxy(int FractalRank,int numbers,double total_mass,vector <double>& masses,double G,
		     vector <double>& posx,vector <double>& posy,vector <double>& posz,
		     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    double delta=0.4;
    double x_off=-10.0+((double)(FractalRank % 4))*0.1;
    double y_off=-5.0+((double)((FractalRank/4) % 4))*0.1;
    double z_off=-2.5+((double)(FractalRank/16))*0.1;
    int number1=pow(((double)numbers)+0.1,1.0/3.0);
    int number2=number1*number1;
    for(int ni=0;ni<numbers;ni++)
      {
	int nx=ni % number1;
	int ny=(ni/number1) % number1;
	int nz=ni/number2;
	posx[ni]=(double)nx*delta+x_off;
	posy[ni]=(double)ny*delta+y_off;
	posz[ni]=(double)nz*delta+z_off;
	velx[ni]=0.0;
	vely[ni]=0.0;
	velz[ni]=0.0;
      }
  }
}
