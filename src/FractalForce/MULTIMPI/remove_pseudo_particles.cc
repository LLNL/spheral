#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void remove_pseudo_particles(Fractal_Memory& mem,Fractal& frac)
  {
    if(!mem.periodic)
      return;
    ofstream& FF=frac.p_file->DUMPS;
    //    ofstream& FF=frac.p_file->FileFractal;
    FF << " entered remove_particles " << frac.get_number_particles() << "\n";
    FF << " Total number of particles entering Gather " << Particle::number_particles << "\n";
    int total=frac.pseudo_particle_list.size();
    for(int particle=0;particle < total;particle++)
      delete frac.pseudo_particle_list[particle];
    frac.pseudo_particle_list.clear();
    FF << " Total number of particles exiting Gather " << Particle::number_particles << "\n";
    FF << " leaving remove_particles " << frac.get_number_particles() << "\n";
  }
}  
