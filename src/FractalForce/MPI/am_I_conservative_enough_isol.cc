#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void am_I_conservative_enough_isol(Fractal_Memory* PFM,vector <double>& masses,double G,
				     vector <double>& xmin,vector <double>& xmax,double correction,
				     vector <double>& posx,vector <double>& posy,vector <double>& posz,
				     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    ofstream& FileEnergy=PFM->p_file->FileEnergy;
    ofstream& FileMom=PFM->p_file->FileMom;
    ofstream& FileP=PFM->p_file->FileParticle;
    int stride=100;
    int NP=PFM->number_particles;
    double dt5=correction*PFM->step_length;
    vector <double>pot(stride);
    vector <double>fx(stride);
    vector <double>fy(stride);
    vector <double>fz(stride);
    double vx,vy,vz;
    double pe=0.0;
    double ke=0.0;
    double p0=0.0;
    double p1=0.0;
    double p2=0.0;
    double m0=0.0;
    double m1=0.0;
    double m2=0.0;
    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	get_field(PFM,ni,many,G,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    if(!I_am_a_real_particle(PFM,nip))
	       continue;
	    vx=velx[nip]+fx[p]*dt5;
	    vy=vely[nip]+fy[p]*dt5;
	    vz=velz[nip]+fz[p]*dt5;
	    pe+=masses[nip]*pot[p];
	    ke+=masses[nip]*(vx*vx+vy*vy+vz*vz);
	    p0+=masses[nip]*vx;
	    p1+=masses[nip]*vy;
	    p2+=masses[nip]*vz;
	    m0+=masses[nip]*(posy[nip]*vz-posz[nip]*vy);
	    m1+=masses[nip]*(posz[nip]*vx-posx[nip]*vz);
	    m2+=masses[nip]*(posx[nip]*vy-posy[nip]*vx);
	  }
      }
    pe*=0.5;
    ke*=0.5;
    double te=pe+ke;
    vector <double>sums(9);
    sums[0]=te;
    sums[1]=pe;
    sums[2]=ke;
    sums[3]=p0;
    sums[4]=p1;
    sums[5]=p2;
    sums[6]=m0;
    sums[7]=m1;
    sums[8]=m2;
    PFM->p_mess->Find_Sum_DOUBLE(sums,9);
    FileEnergy << scientific << PFM->time <<  "\t " << PFM->steps << "\t " << 
      te << "\t " << pe << "\t " << ke << "\t " << sums[0] << "\t" << sums[1] << "\t" << sums[2] << endl;
    FileMom << scientific << PFM->time << "\t " << PFM->steps << "\t " << 
      p0 << "\t " << p1 << "\t " << p2 << "\t " <<
      m0 << "\t " << m1 << "\t " << m2 << "\t " << 
      sums[3] << "\t" << sums[4] << "\t" << sums[5] << "\t" <<
      sums[6] << "\t" << sums[7] << "\t" << sums[8];
  }
}
