#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  typedef deque<double>::iterator _ITD__;
  template <class ForwardIterator>
  void take_a_leap_isol(Fractal_Memory* PFM,ForwardIterator massesb,double G,
			vector <double>& xmin,vector <double>& xmax,
			ForwardIterator posxb,ForwardIterator posyb,
			ForwardIterator poszb,ForwardIterator velxb,
			ForwardIterator velyb,ForwardIterator velzb)
  {
    int stride=100;
    int NP=PFM->number_particles;
    double dt=PFM->step_length;
    fprintf(PFM->p_file->PFFractalMemory," leap isol %10.2E %10.2E %d \n ",PFM->time,dt,NP);
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
	    if(I_am_a_real_particle(PFM,ni+p))
	      {
		*velxb+=fx[p]*dt;
		*velyb+=fy[p]*dt;
		*velzb+=fz[p]*dt;
	      }
	    *posxb+=*velxb*dt;
	    *posyb+=*velyb*dt;
	    *poszb+=*velzb*dt;
	    posxb++;
	    posyb++;
	    poszb++;
	    velxb++;
	    velyb++;
	    velzb++;
	  }
      }
  }
}
namespace FractalSpace
{
  template
  void take_a_leap_isol(Fractal_Memory* PFM,_ITD__ massesb,double G,
			vector <double>& xmin,vector <double>& xmax,
			_ITD__ posxb,_ITD__ posyb,
			_ITD__ poszb,_ITD__ velxb,
			_ITD__ velyb,_ITD__ velzb);
}
namespace FractalSpace
{
  void take_a_leap_isol(Fractal_Memory* PFM,vector <double>& masses,double G,
			vector <double>& xmin,vector <double>& xmax,
			vector <double>& posx,vector <double>& posy,vector <double>& posz,
			vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    //
    int stride=100;
    int NP=PFM->number_particles;
    double dt=PFM->step_length;
    fprintf(PFM->p_file->PFFractalMemory," leap isol %10.2E %10.2E %d \n ",PFM->time,dt,NP);
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
