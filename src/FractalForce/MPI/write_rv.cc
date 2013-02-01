#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void write_rv(const int& step,Fractal& fractal)
  {
    ofstream& FilePos=fractal.p_file->FilePos;
    int nphase=-1;
    int nfield=-1;
    int nparts=fractal.get_number_particles();
    for(int n=0;n<nparts;++n)
      {
	Particle& particle=*fractal.particle_list[n];
	particle.get_phase_field_sizes(nphase,nfield);
	vector <double> pos(3);
	vector <double> vel(3);
	vector <double> pf(4);
	particle.get_pos(pos);
	FilePos << "stepout " << step << "\t" << n << "\t" << fixed << pos[0] << "\t" <<pos[1] << "\t" << pos[2];
	if(nphase >= 6)
	  particle.get_vel(vel);
	else
	  vel.assign(3,-1234.5);
	FilePos << scientific << "\t" << vel[0] << "\t" << vel[1] << "\t" << vel[2];
	FilePos << "\t" << particle.get_mass() ;
	if(nfield >= 4)
	  particle.get_field_pf(pf);
	else
	  pf.assign(4,-1234.5);
	FilePos << scientific << "\t" << pf[0] << "\t" << pf[1] << "\t" << pf[2] << "\t" << pf[3] ;
	FilePos  << endl;
      }
  }
}
