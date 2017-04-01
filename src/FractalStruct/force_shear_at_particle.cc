#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_shear_at_particle(Fractal_Memory& fractal_memory,Fractal& fractal)
  { 
    ofstream& FileFractal=fractal_memory.p_fractal->p_file->FileFractal;
    vector <double> dens(8);
    vector <double> weights(8);
    vector <double> f_xx(8);
    vector <double> f_xy(8);
    vector <double> f_xz(8);
    vector <double> f_yy(8);
    vector <double> f_yz(8);
    vector <double> f_zz(8);
    vector <double> shear(6);
    vector <double> shear2(6,0.0);
    const  double scale=(double)(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
    double fourpiinv=1.0/(16.0*atan(1.0));
    for(int level=0;level <= fractal.get_level_max();level++)
      {
	for(auto &p_group : fractal_memory.all_groups[level])
	  {
	    double d_inv=pow(2.0,p_group->get_level()-fractal.get_level_max());
	    for(auto &p_point : p_group->list_points)
	      {
		Point& point=*p_point;
		bool not_yet=true;
		for(auto &p_part : point.list_particles)
		  {
		    Particle& particle=*p_part;
		    if(!particle.get_real_particle())
		      continue;
		    if(p_group != particle.get_p_highest_level_group())
		      continue;
		    if(not_yet)
		      {
			point.get_field_shear_values(f_xx,f_xy,f_xz,f_yy,f_yz,f_zz);
			not_yet=false;
		      }
		    double d_x,d_y,d_z;
		    vector <double> pos(3);
		    particle.get_pos(pos);
		    point.get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
		    Misc::set_weights(weights,d_x,d_y,d_z);
		    Misc::sum_prod<double>(0,7,1,shear,weights,f_xx,f_xy,f_xz,f_yy,f_yz,f_zz);
		    if(fractal_memory.calc_density_particle)
		      particle.set_density(-(shear[0]+shear[3]+shear[5])*fourpiinv);
		    if(!fractal_memory.calc_shear)
		      continue;
		    for(int i=0;i<6;i++)
		      shear2[i]+=shear[i]*shear[i];
		    double r=0.0;
		    max_predict(fractal_memory,fractal,shear,r);
		    if(fractal_memory.start_up) 
		      particle.set_rad_max(r);
		  }
	      }
	  }
      }
    for(int level=0;level <= fractal.get_level_max();level++)
      for(auto &p_group : fractal_memory.all_groups[level])
	for(auto &p_point : p_group->list_points)
	  p_point->force_shear_point_zero();
    if(fractal_memory.calc_shear)
      {
	double spam=(double)fractal.get_number_particles();
	FileFractal << " shear2 " ;
	for(int ni=0;ni<6;ni++)
	  FileFractal << shear2[ni]/spam << " " ;
	FileFractal << "\n";
      }
  }
}
