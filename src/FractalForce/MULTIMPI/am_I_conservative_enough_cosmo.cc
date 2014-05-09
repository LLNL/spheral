#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void am_I_conservative_enough_cosmo(Fractal_Memory* PFM,vector <double>& masses,double G,
				     vector <double>& xmin,vector <double>& xmax,double correction,
				     vector <double>& posx,vector <double>& posy,vector <double>& posz,
				     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    ofstream& FileEnergy=PFM->p_file->FileEnergy;
    ofstream& FileMom=PFM->p_file->DUMPS;
    //    ofstream& FileMom=PFM->p_file->FileMom;
    //    ofstream& FileP=PFM->p_file->FileParticle;
    int stride=100;
    int NP=PFM->number_particles;
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

    double parad=pow(PFM->arad,PFM->pexp);
    double dpda=PFM->pexp*parad/PFM->arad;
    double parad_half=parad+0.5*PFM->step_length;
    double arad_half=pow(parad_half,1.0/PFM->pexp);
    //    double dadt=arad_half*Fractal_Memory::hubble(arad_half,PFM->omega_start,PFM->lambda_start);
    double v_consta=1.0-2.0*correction*PFM->step_length/PFM->arad/dpda;
    double f_const=PFM->step_length/pow(PFM->arad,4)/Fractal_Memory::hubble(arad_half,PFM->omega_start,PFM->lambda_start)/dpda;
    double v_constb=correction*f_const*(1.0+2.0*PFM->step_length/PFM->arad/dpda);

    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	get_field(PFM,ni,many,G,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    if(!I_am_a_real_particle(PFM,nip))
	       continue;
	    vx=velx[nip]*v_consta+fx[p]*v_constb;
	    vy=vely[nip]*v_consta+fy[p]*v_constb;
	    vz=velz[nip]*v_consta+fz[p]*v_constb;
	    pe+=masses[nip]*pot[p];
	    ke+=masses[nip]*(vx*vx+vy*vy+vz*vz);
	    p0+=masses[nip]*vx;
	    p1+=masses[nip]*vy;
	    p2+=masses[nip]*vz;
	  }
      }
    pe*=0.5;
    ke*=0.5;
    double te=pe+ke;
    vector <double>sums(6);
    sums[0]=te;
    sums[1]=pe;
    sums[2]=ke;
    sums[3]=p0;
    sums[4]=p1;
    sums[5]=p2;
    PFM->p_mess->Find_Sum_DOUBLE(sums,6);
    FileEnergy << scientific << PFM->time <<  "\t " << PFM->steps << "\t " << 
      te << "\t " << pe << "\t " << ke << "\t " << sums[0] << "\t" << sums[1] << "\t" << sums[2] << "\n";
    FileMom << scientific << PFM->time << "\t " << PFM->steps << "\t " << 
      p0 << "\t " << p1 << "\t " << p2 << "\t " <<
      sums[3] << "\t" << sums[4] << "\t" << sums[5] << "\t" << "\n";
  }
}
