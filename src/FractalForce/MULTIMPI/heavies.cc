#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void heavies(Fractal& fractal,Fractal& fractal_ghost)
  {
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    int n_ghost=0;
    fractal_ghost.set_number_particles(n_ghost);
    if(fractal.get_force_max() <= 0.0) return;
    vector <double> spacing(fractal.get_level_max()+1);
    vector <double> offset(fractal.get_level_max()+1);
    vector <double> smaller(fractal.get_level_max()+1);
    for(int lev=0;lev <= fractal.get_level_max();++lev)
      {
	spacing[lev]=pow(2.0,-lev)/((double)fractal.get_grid_length());
	offset[lev]=-0.5*spacing[lev];
	smaller[lev]=pow(8.0,-lev);
      }
    const double log4=log(4.0);
    const double tmp_mass=fractal.get_force_max()/(double)Misc::pow(fractal.get_grid_length(),2);
    const double very_small=1.0e-30;
    for(int p=0;p< fractal.get_number_particles();++p)
      {
	Particle& particle=*fractal.particle_list[p];
	Group* p_group=particle.get_p_highest_level_group();
	if(p_group ==0) continue;
	int lev=p_group->get_level();
	int h_p=(int)floor(log(tmp_mass/(particle.get_mass()+very_small))/log4);
	if(h_p < lev)
	  {
	    n_ghost+=Misc::pow(8,lev-h_p)+1;
	    FileFractal << "n_ghost " << n_ghost << " " << p << "\n";
	  }
      }
    //    assert(n_ghost < fractal.get_number_particles());
    fractal_ghost.set_number_particles(n_ghost);
    fractal_ghost.particle_list.resize(n_ghost);
    if(n_ghost == 0) return;

    //
    Particle* ph;
    try
      {
	ph=new Particle[n_ghost];
      }
    catch(bad_alloc& ba)
      {
	cerr << " Bad Ghost " << n_ghost << " " << ba.what() << endl;
	exit(0);
      }
    int n_h=0; 
    for(int p=0;p< fractal.get_number_particles();++p)
      {
	Particle& particle=*fractal.particle_list[p];
	Group* p_group=particle.get_p_highest_level_group();
	if(p_group ==0) continue;
	int lev=p_group->get_level();
	int h_p=(int)(log(tmp_mass/(particle.get_mass()+very_small))/log4);
	if(h_p < lev) 
	  {
	    const int n_p=Misc::pow(2,lev-h_p);
	    const double delta=spacing[h_p]/(double)n_p;
	    const double off=-delta*0.5*(double)(n_p-1);
	    const double small=smaller[lev-h_p];
	    vector <double> pos(3);
	    particle.get_pos(pos);
	    double mm=particle.get_mass()*small;
	    for(int p_z=0;p_z < n_p;++p_z)
	      {
		for(int p_y=0;p_y < n_p;++p_y)
		  {
		    for(int p_x=0;p_x < n_p;++p_x)
		      {
			fractal_ghost.particle_list[n_h]=&ph[n_h];
			Particle& p_ghost=ph[n_h];
			p_ghost.set_pos(pos[0]+(double)p_x*delta+off,pos[1]+(double)p_y*delta+off,pos[2]+(double)p_z*delta+off);
			p_ghost.set_mass(mm);
			++n_h;
		      }
		  }
	      }
	    fractal_ghost.particle_list[n_h]=&ph[n_h];
	    Particle& p_ghost=ph[n_h];
	    p_ghost.set_pos(pos);
	    p_ghost.set_mass(-particle.get_mass());
	    ++n_h;
	  }
      }
    assert(n_h==n_ghost);
  }
}
