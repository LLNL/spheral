#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_particle_highest(Group& group, Fractal& fractal)
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
    for(list<Point*>::const_iterator point_itr=group.list_points_highest.begin();point_itr !=group.list_points_highest.end();++point_itr)
      {
	Point& point=**point_itr;
	bool not_yet=true;
	//
	for(list<Particle*>::const_iterator particle_itr=point.list_particles_highest.begin();particle_itr !=point.list_particles_highest.end();++particle_itr)
	  {
	    Particle& particle=**particle_itr;
	    if(not_yet)
	      {
		point.get_field_values(pott,f_x,f_y,f_z);
		not_yet=false;
	      }
	    particle.get_pos(pos);
	    point.get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
	    Misc::set_weights(weights,d_x,d_y,d_z);
	    Misc::sum_prod<double>(0,7,1,sum_pf,weights,pott,f_x,f_y,f_z);
	    particle.set_field_pf(sum_pf);
	    if(sum_pf[0]*sum_pf[1]*sum_pf[2]*sum_pf[3] ==0.0)
	      {
		cout << "part forces b " << &particle;
		vector <double> pos(3);
		vector <double> field(4);
		particle.get_pos(pos);
		particle.get_field_pf(field);
		cout << " " << pos[0];
		cout << " " << pos[1];
		cout << " " << pos[2];
		cout << " " << field[0];
		cout << " " << field[1];
		cout << " " << field[2];
		cout << " " << field[3];
		cout << endl;
		for(int ii=0;ii < 8;ii++)
		  cout << pott[ii] << " " << f_x[ii] << " " << f_y[ii] << " " << f_z[ii] << endl;
	      }
	  }
      }
  }
}
