#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_particle_sharp(Group& group, Fractal& fractal)
  { 
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    //    FileFractal << " enter force at particle sharp" << "\n";
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
    for(auto p : group.list_points)
      {
	bool not_yet=true;
	//
	for(auto part : p->list_particles)
	  {
	    if(!part->get_real_particle())
	      continue;
	    if(part->get_p_highest_level_group() != 0)
	      {
		if(p_group == part->get_p_highest_level_group())
		  {
		    if(not_yet)
		      {
			p->get_field_values(pott);
			not_yet=false;
		      }
		    //		    FileFractal << "sharp " << &group << " " << &point << " " << &particle <<"\n";
		    part->get_pos(pos);
		    p->get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
		    Misc::set_weights(weights_p,weights_x,weights_y,weights_z,d_x,d_y,d_z);
		    Misc::sum_prod_p_sharp(0,7,1,sum_pf,weights_p,weights_x,weights_y,weights_z,pott);
		    part->set_field_pf(sum_pf,scale_force);
		    if(sum_pf[0]*sum_pf[1]*sum_pf[2]*sum_pf[3] ==0.0)
		      part->dump(FileFractal,pott,f_x,f_y,f_z);
		    //		    part->dump(FileFractal);
		  }
	      }
	    else
	      {
		part->dump(FileFractal);
		part->set_field_pf(0.0);
	      }
	  }
      }
    //    FileFractal << " exit force at particle sharp" << "\n";
  }
}
