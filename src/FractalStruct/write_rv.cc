#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void write_rv(const int& step,Fractal& fractal)
  {
    //    ofstream& FilePos=fractal.p_file->FilePos;
    FILE* PFPos=fractal.p_file->PFPos;
    int nphase=-1;
    int nfield=-1;
    int nparts=fractal.get_number_particles();
    double base_mass=fractal.get_base_mass();
    bool period=fractal.get_periodic();
    int mlevel=-1;
    double amlevel=0.0;
    for(int n=0;n<nparts;++n)
      {
	Particle& particle=*fractal.particle_list[n];
	particle.get_phase_field_sizes(nphase,nfield);
	int plevel=particle.get_highest_level();
	vector <double> pos(3);
	vector <double> vel(3);
	vector <double> pf(4);
	particle.get_pos(pos);
	fprintf(PFPos," Out %d %7d %10.6E %10.6E %10.6E",step,n,pos[0],pos[1],pos[2]);
	if(nphase >= 6)
	  particle.get_vel(vel);
	else
	  vel.assign(3,-1234.5);
	fprintf(PFPos," %10.6E %10.6E %10.6E",vel[0],vel[1],vel[2]);
	double m=particle.get_mass();
	//	FilePos << "\t" << m ;
	if(period)
	  {
	    if(m > 0.0)
	      {
		amlevel=base_mass/m;
		mlevel=(amlevel+0.01);
	      }
	    else
	      mlevel=-1;
	    fprintf(PFPos," M %d",mlevel);
	  }
	fprintf(PFPos," L %d",plevel);
	if(nfield >= 4)
	  particle.get_field_pf(pf);
	else
	  pf.assign(4,-1234.5);
	fprintf(PFPos," %10.6E %10.6E %10.6E %10.6E \n",pf[0],pf[1],pf[2],pf[3]);
      }
    //    fflush(PFPos);
  }
}
