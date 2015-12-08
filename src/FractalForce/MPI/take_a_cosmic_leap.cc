#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void take_a_cosmic_leap(Fractal_Memory* PFM,vector <double>& masses,double G,
			vector <double>& xmin,vector <double>& xmax,
			vector <double>& posx,vector <double>& posy,vector <double>& posz,
			vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    int stride=100;
    int NP=PFM->number_particles;
    vector <double>pot(stride);
    vector <double>fx(stride);
    vector <double>fy(stride);
    vector <double>fz(stride);

    double parad=pow(PFM->arad,PFM->pexp);
    double dpda=PFM->pexp*parad/PFM->arad;
    double parad_half=parad+0.5*PFM->step_length;
    double arad_half=pow(parad_half,1.0/PFM->pexp);
    double dadt=arad_half*Fractal_Memory::hubble(arad_half,PFM->omega_start,PFM->lambda_start);
    double v_const=1.0-2.0*PFM->step_length/PFM->arad/dpda;
    double f_const=PFM->step_length/pow(PFM->arad,4)/Fractal_Memory::hubble(arad_half,PFM->omega_start,PFM->lambda_start)/dpda;
    double dt=PFM->step_length/dadt/dpda;

    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	get_field(PFM,ni,many,G,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    if(!I_am_a_real_particle(PFM,nip))
	      continue;
	    velx[nip]=velx[ni]*v_const+fx[p]*f_const;
	    vely[nip]=vely[ni]*v_const+fy[p]*f_const;
	    velz[nip]=velz[ni]*v_const+fz[p]*f_const;
	    posx[nip]+=velx[nip]*dt;
	    posy[nip]+=vely[nip]*dt;
	    posz[nip]+=velz[nip]*dt;
	  }
      }
    PFM->time+=PFM->step_length/dadt/dpda;
    if(!PFM->start_up)
      parad+=PFM->step_length;
    PFM->arad=pow(parad,1.0/PFM->pexp);
  }
}
