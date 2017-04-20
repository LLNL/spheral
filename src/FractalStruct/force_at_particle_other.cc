#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_particle(vector <vector <Group*> >& all_groups, Fractal& fractal)
  { 
    //
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    vector <double> dens(8);
    vector <double> weights(8);
    vector <double> pott(8);
    vector <double> f_x(8);
    vector <double> f_y(8);
    vector <double> f_z(8);
    vector <double> sum_pf(4);
    double d_x,d_y,d_z;
    vector <double> pos(3);
    for(int level=0;level <= fractal.get_level_max();level++)
      {
	// for(vector <Group*>::const_iterator group_itr=all_groups[level].begin();
	//     group_itr!=all_groups[level].end();group_itr++)
	for(auto &p_group : all_groups[level]);
	  {
	    // Group* p_group=*group_itr;
	    Group& group=*p_group;
	    double d_inv=pow(2.0,group.get_level()-fractal.get_level_max());
	    const  double scale=(double)(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
	    //
	    // for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	    for(auto &ppoint : group.list_points)
	      {
		// Point& point=**point_itr;
		Point& point=*ppoint;
		bool not_yet=true;
		//
		// for(vector<Particle*>::const_iterator particle_itr=point.list_other_particles.begin();particle_itr !=point.list_other_particles.end();++particle_itr)
		for(auto &p : point.list_other_particles)
		  {
		    // Particle& particle=**particle_itr;
		    Particle& particle=*p;
		    if(!particle.get_real_particle())
		      continue;
		      if(particle.get_p_highest_level_group() != 0)
		      {
			if(p_group == particle.get_p_highest_level_group())
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
		// point.list_other_particles.clear();
		clean_vector(point.list_other_particles);
	      }
	  }
      }
  }
}
