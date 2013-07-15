#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_particle_sharp(Group& group, Fractal& fractal)
  { 
    ofstream& FileFractal=fractal.p_file->FileFractal;
    FileFractal << " enter force at particle sharp" << endl;
    vector <double> weights_p(8);
    vector <double> weights_x(8);
    vector <double> weights_y(8);
    vector <double> weights_z(8);
    vector <double> pott(8);
    vector <double> f_x(8);
    vector <double> f_y(8);
    vector <double> f_z(8);
    vector <double> sum_pf(4);
    double d_x,d_y,d_z;
    vector <double> pos(3);
    Group* p_group=&group;
    double d_inv=pow(2.0,group.get_level()-fractal.get_level_max());
    const  double scale=(double)(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
    double scale_force=(double)(fractal.get_grid_length()*Misc::pow(2,group.get_level()));
    //
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	bool not_yet=true;
	//
	for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
	  {
	    Particle& particle=**particle_itr;
	    if(!particle.get_real_particle())
	      continue;
	    if(particle.get_p_highest_level_group() != 0)
	      {
		if(p_group == particle.get_p_highest_level_group())
		  {
		    if(not_yet)
		      {
			point.get_field_values(pott);
			not_yet=false;
		      }
		    //		    FileFractal << "sharp " << &group << " " << &point << " " << &particle <<endl;
		    particle.get_pos(pos);
		    point.get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
		    Misc::set_weights(weights_p,weights_x,weights_y,weights_z,d_x,d_y,d_z);
		    Misc::sum_prod_p_sharp(0,7,1,sum_pf,weights_p,weights_x,weights_y,weights_z,pott);
		    particle.set_field_pf(sum_pf,scale_force);
		    if(sum_pf[0]*sum_pf[1]*sum_pf[2]*sum_pf[3] ==0.0)
		      particle.dump(FileFractal,pott,f_x,f_y,f_z);
		    //		    particle.dump(FileFractal);
		  }
	      }
	    else
	      {
		particle.dump(FileFractal);
		particle.set_field_pf(0.0);
	      }
	  }
      }
    FileFractal << " exit force at particle sharp" << endl;
  }
}
