#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void buffer_points(Group& group, Fractal& fractal)
  {
    fractal.timing(-1,10);
    ofstream& FileFractal=fractal.p_file->DUMPS;
    int padd=fractal.get_padding();
    double scale=static_cast<double>(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
    int width=Misc::pow(2,fractal.get_level_max()-group.get_level()-1);
    if(padd > 0)
      {
	for(auto &p : group.p_list_really_high)
	  {
	    Point* p_point=p;
	    for(int p=0;p < padd;p++)
	      for(int ni=0;ni<3;ni++)
		p_point=p_point->get_point_ud_0(ni*2,0);
	    Point* p_point_z=p_point;
	    for(int pad_z=-padd;pad_z<=padd;++pad_z)
	      {
		Point* p_point_y=p_point_z;
		for(int pad_y=-padd;pad_y <= padd;++pad_y)
		  {
		    Point* p_point_x=p_point_y; 
		    for(int pad_x=-padd;pad_x <= padd;++pad_x)
		      {
			p_point_x->set_it_is_high(true);
			if(pad_x < padd)
			  p_point_x=p_point_x->get_point_up_x_0(); 
		      }
		    if(pad_y < padd)
		      p_point_y=p_point_y->get_point_up_y_0();
		  }
		if(pad_z < padd)
		  p_point_z=p_point_z->get_point_up_z_0();
	      }
	  }
      }
    else if( padd==-1)
      {
	unsigned int min_number=fractal.get_minimum_number();
	vector <unsigned int> numbers(8,0);
	vector <double> pos(3);
	int p_x,p_y,p_z,n_x,n_y,n_z;
	for(auto &p : group.p_list_really_high)
	  {
	    Point& point=*p;
	    if(point.list_particles.size() < min_number) continue;
	    numbers.assign(8,0);
	    for(auto &part : point.list_particles)
	      {
	    	Particle& particle=*part;
		particle.get_pos(pos);
		point.get_pos_point(p_x,p_y,p_z);
		n_x=(static_cast<int>(floor(pos[0]*scale))-p_x)/width;
		n_y=(static_cast<int>(floor(pos[1]*scale))-p_y)/width;
		n_z=(static_cast<int>(floor(pos[2]*scale))-p_z)/width;
		if(n_x > 1 || n_y > 1 || n_z > 1 || n_x < 0 || n_y < 0 || n_z < 0)
		  {
		    fractal.p_mess->Full_Stop();
		    FileFractal << " buffer badda " << endl;
		    point.dump();
		    particle.dump(*point.p_FILE);
		    FileFractal << " buffer baddb " << endl;
		  }
		assert(n_x==0 || n_x==1);
		assert(n_y==0 || n_y==1);
		assert(n_z==0 || n_z==1);
		assert(pos[0]*scale >= p_x);
		assert(pos[1]*scale >= p_y);
		assert(pos[2]*scale >= p_z);
		int corner=n_x+2*n_y+4*n_z;
		numbers[corner]++;
	      }
	    for(int corner=0;corner < 8;++corner)
	      {
		if(numbers[corner] >= min_number)
		  list_buffer(point,corner);
	      }
	  }
      }
    if(!fractal.get_periodic())
      {
	int zoom=fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max());
	vector <int>pose(3);
	for(auto &p : group.list_points)
	  {
	    if(!p->get_it_is_high())
	      continue;
	    p->get_pos_point(pose);
	    if(
	       pose[0] == 0 || pose[1] == 0 || pose[2] == 0 ||
	       pose[0] == zoom || pose[1] == zoom || pose[2] == zoom
	       )
	      p->set_it_is_high(false);
	  }
      }
    //--------------------------------------------------------------------------------------------------------------------------------
    // make sure edge has no high points for isolated BC
    //--------------------------------------------------------------------------------------------------------------------------------
    int n_h=0;
    for(auto &p : group.list_points)
      {
	p->set_passive_low();
	if(p->get_it_is_high()) n_h++;
      }
    group.set_number_high_points(n_h);
    group.p_list_really_high.clear();
    fractal.timing(1,10);
  }
}
