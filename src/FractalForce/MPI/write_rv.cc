#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void write_rv(const int& step,Fractal& fractal)
  {
    static ofstream FilePos;
    if(!FilePos.is_open())
      FilePos.open("pos.d");
    //      FilePos.open("/p/lscratcha/jensv/pos.d");
    FilePos.precision(7);
    for(int n=0;n<fractal.get_number_particles();++n)
      {
	Particle& particle=*fractal.particle_list[n];
	vector <double> pos(3);
	vector <double> vel(3);
	particle.get_phase(pos,vel);
	FilePos << "stepout " << step << "\t" << n << "\t" << fixed << pos[0] << "\t" <<pos[1] << "\t" << pos[2];
	FilePos << scientific << "\t" << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t" << particle.get_mass() ;
	FilePos << "\t" << particle.get_density();
	FilePos  << endl;
      }
  }
}
