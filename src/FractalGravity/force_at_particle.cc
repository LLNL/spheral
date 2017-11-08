#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_particle(Group& group, Fractal& fractal,const bool& conserve)
  { 
    //
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
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(point.list_particles.empty()) continue;
	bool not_yet=true;
	//
	for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
	  {
	    Particle& particle=**particle_itr;
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
		      particle.dump(pott,f_x,f_y,f_z);
		  }
	      }
	    else
	      {
		particle.dump();
		particle.set_field_pf(0.0);
	      }
	  }
      }
  }
}
