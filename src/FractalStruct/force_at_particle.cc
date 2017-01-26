#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_particle(Group& group, Fractal& fractal,const bool& conserve)
  { 
    fractal.timing(-1,8);
    ofstream& FileFractal=fractal.p_file->DUMPS;
    vector <double> dens(8);
    vector <double> weights(8);
    vector <double> pott(8);
    vector <double> f_x(8);
    vector <double> f_y(8);
    vector <double> f_z(8);
    vector <double> sum_pf(4);
    Group* p_group=&group;
    double d_x,d_y,d_z;
    vector <double> pos(3);
    double d_inv=pow(2.0,group.get_level()-fractal.get_level_max());
    const  double scale=(double)(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
    //
    for(auto &p : group.list_points)
      {
	Point& point=*p;
	if(point.list_particles.empty()) continue;
	bool not_yet=true;
	for(auto &part : point.list_particles)
	  {
	    Particle& particle=*part;
	    if(!particle.get_real_particle())
	      continue;
	    if(particle.get_p_highest_level_group() != 0)
	      {
		if(conserve || p_group == particle.get_p_highest_level_group())
		  {
		    if(not_yet)
		      {
			point.get_field_values(pott,f_x,f_y,f_z);
			not_yet=false;
		      }
		    //
		    particle.get_pos(pos);
		    point.get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
		    Misc::set_weights(weights,d_x,d_y,d_z);
		    Misc::sum_prod<double>(0,7,1,sum_pf,weights,pott,f_x,f_y,f_z);
		    particle.set_field_pf(sum_pf);
		    if(sum_pf[0]*sum_pf[1]*sum_pf[2]*sum_pf[3] ==0.0)
		      particle.dump(FileFractal,pott,f_x,f_y,f_z);
		  }
	      }
	    else
	      {
		particle.dump(FileFractal);
		particle.set_field_pf(0.0);
	      }
	  }
      }
    fractal.timing(1,8);
  }
}
