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
    // if(Mess::IAMROOT)
    //   {
    // 	cerr << " entered remove_particles " << frac.get_number_particles() << "\n";
    // 	cerr << " Total number of particles entering Gather " << Particle::number_particles << "\n";
    //   }
    FF << " entered remove_particles " << frac.get_number_particles() << "\n";
    FF << " Total number of particles entering Gather " << Particle::number_particles << "\n";
    for(auto &p : frac.pseudo_particle_list)
      delete p;
    // frac.particle_list.resize(frac.particle_list.size()-frac.pseudo_particle_list.size());
    // frac.set_number_particles(frac.particle_list.size());
    clean_deque(frac.pseudo_particle_list);
    // if(Mess::IAMROOT)
    //   {
    // 	cerr << " Total number of particles exiting Gather " << Particle::number_particles << "\n";
    // 	cerr << " leaving remove_particles " << frac.get_number_particles() << "\n";
    //   }
    FF << " Total number of particles exiting Gather " << Particle::number_particles << "\n";
    FF << " leaving remove_particles " << frac.get_number_particles() << "\n";
  }
}  
