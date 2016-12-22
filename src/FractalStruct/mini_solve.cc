#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void mini_solve(Fractal_Memory& mem,Group* pg)
  {
    static const double pi = 4.0*atan(1.0);
    int lps=pg->list_points.size();
    if(lps > mem.hypre_max_node_load)
      return;
    double length=mem.p_fractal->get_grid_length();
    double gc=4.0*pi/(length*length)*pow(4.0,-pg->get_level());
    if(lps == 27)
      {
	for(auto &p : pg->list_points)
	  {
	    if(p->get_inside())
	      {
		mini_solve1(p,gc);
		break;
	      }
	  }
      }
    else if(lps == 45)
      {
    	vector<Point*>found;
    	for(auto &p : pg->list_points)
    	  {
    	    if(found.size() == 2)
    	      break;
    	    if(p->get_inside())
    	      found.push_back(p);
    	  }
	assert(found.size() == 2);
	mini_solve2(found[0],found[1],gc);
      }
  }
}
namespace FractalSpace
{
  void mini_solve1(Point* p,const double& gc)
  {
    double sum=0.0;
    for(int ni=0;ni<6;ni++)
      sum+=p->get_point_ud_0(ni)->get_potential_point();
    p->set_potential_point((sum-gc*p->get_density_point())/6.0);
  }
}
namespace FractalSpace
{
  void mini_solve2(Point* pa,Point* pb,const double& gc)
  {
    double suma=-pb->get_potential_point();
    double sumb=-pa->get_potential_point();
    for(int ni=0;ni<6;ni++)
      {
	suma+=pa->get_point_ud_0(ni)->get_potential_point();
	sumb+=pb->get_point_ud_0(ni)->get_potential_point();
      }
    double AB=(suma-gc*pa->get_density_point())/35.0;
    double CD=(sumb-gc*pb->get_density_point())/35.0;
    pa->set_potential_point(6.0*AB+CD);
    pb->set_potential_point(AB+6.0*CD);
  }
}
