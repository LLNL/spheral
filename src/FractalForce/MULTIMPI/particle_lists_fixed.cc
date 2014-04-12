#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void particle_lists_fixed(vector <vector <Group*> >& all_groups,Fractal& fractal,Misc& misc)
  {
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    FileFractal << "enter particle lists " << fractal.get_level_max() << " " << fractal.get_number_particles() << "\n";
    if(fractal.get_level_max() < 1) return;
    if(fractal.get_number_particles() < 1) return;
    const double grid_length=fractal.get_grid_length();
    const double grid_multiplier=(double)misc.grid_multiply;
    vector <double> pos(3);
    Group& group=**all_groups[0].begin();
    for(int particle=0;particle < fractal.get_number_particles();++particle)
      {
	FileFractal << " haha 1 " << "\n";
	Particle& p=*fractal.particle_list[particle];
	p.get_pos(pos);
	int n_x=(int)floor(pos[0]*grid_length);
	int n_y=(int)floor(pos[1]*grid_length);
	int n_z=(int)floor(pos[2]*grid_length);
	int pp=fractal.where_1(n_x,n_y,n_z);
	bool doit=pp >= 0;
	if(doit)
	  {
	    //
	    FileFractal << particle << " " << n_x << " " << n_y << " " << n_z << " " << pp << "\n";
	    //
	    Point& point=*(group.list_points[pp]);
	    //
	    FileFractal << "point " << pp << " " << &point << "\n";
	    FileFractal << point.list_other_particles.size() << "\n";
	    //
	    point.list_other_particles.push_back(&p);
	    p.set_p_highest_level_group(misc.p_group_0);
	    FileFractal << " haha 2 " << "\n";
	  }
	else
	  {
	    FileFractal << "particle outside " << pos[0] << " " << pos[1] << " " << pos[2] << " " << particle << "\n";
	  }
      }
    for(int lev=1;lev <= fractal.get_level_max();++lev)
      {
	const int d_inv=Misc::pow(2,fractal.get_level_max()-lev);
	for(vector <Group*>::const_iterator group_itr=all_groups[lev].begin();
	    group_itr!=all_groups[lev].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point& point=**point_itr;
		if(point.get_real_pointer() == 0)
		  {
		    Point& m_point=*(point.get_point_pointer());
		    for(vector <Particle*>::const_iterator p_itr=m_point.list_other_particles.begin();p_itr!=m_point.list_other_particles.end();++p_itr)
		      {
			Particle& particle=**p_itr;
			vector<double> pos(3);
			particle.get_pos(pos);
			int p_x=((int)floor(pos[0]*grid_multiplier)-point.get_pos_point_x())/d_inv;
			int p_y=((int)floor(pos[1]*grid_multiplier)-point.get_pos_point_y())/d_inv;
			int p_z=((int)floor(pos[2]*grid_multiplier)-point.get_pos_point_z())/d_inv;
			assert(p_x == 0 || p_x == 1);
			assert(p_y == 0 || p_y == 1);
			assert(p_z == 0 || p_z == 1);
			Point* p_p=&point;
			if(p_x == 1) p_p=p_p->get_point_up_x_0();
			if(p_y == 1) p_p=p_p->get_point_up_y_0();
			if(p_z == 1) p_p=p_p->get_point_up_z_0();
			p_p->list_other_particles.push_back(&particle);
		      }
		  }
	      }
	  }
      }
    FileFractal << "leaving particle lists " << "\n";
  }
}
