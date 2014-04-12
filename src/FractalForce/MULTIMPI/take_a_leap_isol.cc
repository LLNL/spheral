#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void take_a_leap_isol(Fractal_Memory* PFM,vector <double>& masses,double G,
			vector <double>& xmin,vector <double>& xmax,
			vector <double>& posx,vector <double>& posy,vector <double>& posz,
			vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    //    ofstream& FFM=PFM->p_file->FileFractalMemory;
    //
    //    double xcen=10.0;
    //    double ycen=-18.0;
    //    double zcen=11.0;
    //    double rmax=20.0;
    //    double slope=-0.5;
    //    double slope3=slope+3.0;
    //    double myG=3.141592;
    //    double mtotal=1.0e9;
    //    double myGm=myG*mtotal;
    //
    int stride=100;
    int NP=PFM->number_particles;
    double dt=PFM->step_length;
    //    FFM << " leap isol " << PFM->time << " " << dt << " " << NP << "\n";
    fprintf(PFM->p_file->PFFractalMemory," leap isol 10.2E 10.2E %d \n ",PFM->time,dt,NP);
    PFM->time+=dt;
    vector <double>pot(stride);
    vector <double>fx(stride);
    vector <double>fy(stride);
    vector <double>fz(stride);
    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	get_field(PFM,ni,many,G,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    if(I_am_a_real_particle(PFM,nip))
	      {
		velx[nip]+=fx[p]*dt;
		vely[nip]+=fy[p]*dt;
		velz[nip]+=fz[p]*dt;
	      }
	    posx[nip]+=velx[nip]*dt;
	    posy[nip]+=vely[nip]*dt;
	    posz[nip]+=velz[nip]*dt;
	  }
      }
  }
}
