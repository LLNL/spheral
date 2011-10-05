#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main()
{
  Fractal* p_fractal=new Fractal;
  Fractal& fractal=*p_fractal;
  //
  cout << "here ", cout << endl;
  fractal.set_number_particles(262144);
  fractal.set_grid_length(64); 
  fractal.set_periodic(false);
  fractal.set_minimum_number(8);
  fractal.set_force_max(1.0e5);
  fractal.set_level_max(5);
  fractal.set_padding(0);
  fractal.set_tweaks(2);
  cout << "parameters " << fractal.get_number_particles() << " "
       << fractal.get_grid_length() << " "
       << fractal.get_periodic() << " "
       << fractal.get_minimum_number() << endl;
  cout << fractal.get_level_max() << " "
       << fractal.get_padding() << " "
       << fractal.get_tweaks() << endl;
  //
  fractal.resize_all();
  double pi=3.141592653589793;
  int seed=12345;
  srand(seed);
  double rand_max=(double)RAND_MAX;
  for(int n=0;n<fractal.get_number_particles();++n)
    {
      double r=0.4*(double)(rand())/rand_max;
      double phi=2.0*pi*(double)(rand())/rand_max;
      double ct=2.0*(double)(rand())/rand_max-1.0;
      double st=sqrt(1.0-ct*ct);
      fractal.set_pos_x(n,r*st*cos(phi)+0.5);
      fractal.set_pos_y(n,r*st*sin(phi)+0.5);
      fractal.set_pos_z(n,r*ct+0.53);
      //      cout << "pos " << n << " " << fractal.get_pos_x(n) << " " << fractal.get_pos_y(n) << " " << fractal.get_pos_z(n) << endl;
    }
  int n=fractal.get_number_particles()-1;
  double pos_x_0=0.5;
  double pos_y_0=0.5;
  double pos_z_0=0.53;
  fractal.set_pos_x(n,pos_x_0);
  fractal.set_pos_y(n,pos_y_0);
  fractal.set_pos_z(n,pos_z_0);
  fractal.set_particle_mass(n,1.0);
  fractal.set_number_masks(0);
  fractal.set_level_mask(0,0,0);
  fractal.set_pos_mask(0,0.5,0.5,0.5,2.0,2.0,2.0);
  fractal.set_level_mask(1,2,0);
  fractal.set_pos_mask(1,0.5,0.5,0.53,0.2,0.2,0.2);
  fractal.set_level_mask(2,5,1);
  fractal.set_pos_mask(2,0.5,0.5,0.53,0.1,0.1,0.1);
  for(int ii =0;ii < 1 ;++ii)
    {
      fractal.timing(-2,0);
      fractal.timing(-1,29);
      fractal_gravity(fractal);
      fractal.timing(1,29);
      fractal.timing(0,0);
      dump_all_particles(fractal);
    }
  delete p_fractal;
  p_fractal=0;
  return 1;
}
void dump_all_particles(Fractal& fractal)
{
  if(1) return;
  /*
  int n0=fractal.get_number_particles()-1;
  double x0=fractal.get_pos_x(n0);
  double y0=fractal.get_pos_y(n0);
  double z0=fractal.get_pos_z(n0);
  */
  for(int n=0;n<fractal.get_number_particles();++n)
    {
      double dx=fractal.get_pos_x(n)-0.5;
      double dy=fractal.get_pos_y(n)-0.5;
      double dz=fractal.get_pos_z(n)-0.53;
      double fx=fractal.get_force_x(n);
      double fy=fractal.get_force_y(n);
      double fz=fractal.get_force_z(n);
      double potty=fractal.get_potential(n);
      /*
      if(dx < -0.5) dx+=1.0;
      if(dy < -0.5) dy+=1.0;
      if(dz < -0.5) dz+=1.0;
      if(dx > 0.5) dx-=1.0;
      if(dy > 0.5) dy-=1.0;
      if(dz > 0.5) dz-=1.0;
      */
      double r=sqrt(dx*dx+dy*dy+dz*dz+1.0e-30);
      double f=sqrt(fx*fx+fy*fy+fz*fz);
      cout << "output " << n << " " <<  dx << " " << dy << " " << dz << " " << r << " " << potty << " " << f << " " << fx << " " << fy << " " << fz << endl;
    }
}
void dump_group(Group& group,Misc& misc)
{
  if(group.get_level() == 0) return;
  for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
    {
      Point& point=**point_itr;
      double p_x=(double)point.get_pos_point_x()/(double)misc.grid_multiply-0.5;
      double p_y=(double)point.get_pos_point_y()/(double)misc.grid_multiply-0.5;
      double p_z=(double)point.get_pos_point_z()/(double)misc.grid_multiply-0.53;
      /*
      if(p_x < -0.5) p_x+=1.0;
      if(p_x > 0.5) p_x-=1.0;
      if(p_y < -0.5) p_y+=1.0;
      if(p_y > 0.5) p_y-=1.0;
      if(p_z < -0.5) p_z+=1.0;
      if(p_z > 0.5) p_z-=1.0;
      */
      double r=sqrt(p_x*p_x+p_y*p_y+p_z*p_z+1.0e-20);
      double fx=point.get_force_point_x();
      double fy=point.get_force_point_y();
      double fz=point.get_force_point_z();
      double f=sqrt(fx*fx+fy*fy+fz*fz);
      double p0=-1.0/r;
      double f0=1.0/r/r;
      double fx0=-p_x/(r*r*r);
      double fy0=-p_y/(r*r*r);
      double fz0=-p_z/(r*r*r);
      cout << "pot_at_point " << " " << &group << " lev " << group.get_level() << " "  << &point << " " << p_x << " " << p_y << " " << p_z << " " << r 
	   << " " << -point.get_potential_point() << " " << -p0 << " " << f << " " << f0 << " " 
	   << fx << " " << fx0 << " " << fy << " " << fy0 << " " << fz << " " << fz0 << "  ins " << point.get_inside() << endl;
    }
}
