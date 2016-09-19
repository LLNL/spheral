#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void remove_pseudo_particles(Fractal_Memory& mem,Fractal& frac)
  {
    if(!mem.periodic)
      return;
    frac.p_file->FileFractal << " entered remove_particles " << frac.get_number_particles() << endl;
    int counterR=0;
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	if(P->get_real_particle())
	  {
	    if(counterR < particle)
	      frac.particle_list[counterR]=P;
	    counterR++;
	  }
	else
	  delete P;
      }
    frac.particle_list.resize(counterR);
    frac.set_number_particles(counterR);
    frac.p_file->FileFractal << " leaving remove_particles " << frac.get_number_particles() << endl;
  }
}  
