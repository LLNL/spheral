#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void particle_lists(vector <vector <Group*> >& all_groups,Fractal& fractal,Fractal& fractal_other,Misc& misc)
  {
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    FileFractal << "enter particle lists " << fractal.get_level_max() << " " << fractal_other.get_number_particles() << "\n";
    if(fractal.get_level_max() < 1) return;
    if(fractal_other.get_number_particles() < 1) return;
    const double grid_length=fractal.get_grid_length();
    const double grid_multiplier=(double)misc.grid_multiply;
    Group& group=*all_groups[0][0];
    vector <double> pos(3);
    for(int particle=0;particle < fractal_other.get_number_particles();++particle)
      {
	Particle* p_part=fractal_other.particle_list[particle];
	p_part->get_pos(pos);
	int n_x=(int)floor(pos[0]*grid_length);
	int n_y=(int)floor(pos[1]*grid_length);
	int n_z=(int)floor(pos[2]*grid_length);
	int p=fractal.where_1(n_x,n_y,n_z);
	assert(p >= 0);
	Point& point=*(group.list_points[p]);
	point.list_particles.push_back(p_part);
	point.list_other_particles.push_back(p_part);
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
		Point* p_point=*point_itr;
		if(p_point->get_real_pointer() == 0)
		  {
		    Point& m_point=*(p_point->get_point_pointer());
		    if(m_point.list_other_particles.empty()) continue;
		    for(vector <Particle*>::const_iterator p_itr=m_point.list_other_particles.begin();p_itr!=m_point.list_other_particles.end();++p_itr)
		      {
			Particle* p_particle=*p_itr;
			p_particle->get_pos(pos);
			int p_x=((int)floor(pos[0]*grid_multiplier)-p_point->get_pos_point_x())/d_inv;
			int p_y=((int)floor(pos[1]*grid_multiplier)-p_point->get_pos_point_y())/d_inv;
			int p_z=((int)floor(pos[2]*grid_multiplier)-p_point->get_pos_point_z())/d_inv;
			assert(p_x == 0 || p_x == 1);
			assert(p_y == 0 || p_y == 1);
			assert(p_z == 0 || p_z == 1);
			Point* p_p=p_point;
			if(p_x == 1) p_p=p_p->get_point_up_x_0();
			if(p_y == 1) p_p=p_p->get_point_up_y_0();
			if(p_z == 1) p_p=p_p->get_point_up_z_0();
			p_p->list_particles.push_back(p_particle);
			p_p->list_other_particles.push_back(p_particle);
			p_particle->set_p_highest_level_group(0);
			p_particle->set_highest_level(-1);
		      }
		  }
	      }
	  }
      }
    for(int lev=0;lev <= fractal.get_level_max();++lev)
      {
	for(vector <Group*>::const_iterator group_itr=all_groups[lev].begin();
	    group_itr!=all_groups[lev].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point* p_point=*point_itr;
		p_point->list_other_particles.clear();
	      }
	  }
      }

    FileFractal << "leaving particle lists " << "\n";
  }
}
