#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool mini_solve(Fractal_Memory& mem,Group* pg)
  {
    static const double pi = 4.0*atan(1.0);
    bool done=false;
    int lps=pg->list_points.size();
    if(lps > mem.hypre_min_node_load)
      return false;
    double length=mem.p_fractal->get_grid_length();
    double gc=4.0*pi/(length*length)*pow(4.0,-pg->get_level());
    if(lps == 27)
      {
	int count=0;
	for(auto &p : pg->list_points)
	  if(p->get_inside())
	    {
	      mini_solve1(p,gc);
	      count++;
	    }
	done=count==1;
      }
    else if(lps == 45)
      {
    	vector<Point*>found;
    	for(auto &p : pg->list_points)
	  if(p->get_inside())
	    found.push_back(p);
	done=found.size() == 2;
	if(done)
	  mini_solve2(found[0],found[1],gc);
      }
    else if(lps == 63)
      {
    	vector<Point*>found;
    	for(auto &p : pg->list_points)
	  if(p->get_inside())
	    found.push_back(p);
	done=found.size() == 3;
	if(done)
	  mini_solve3(found,gc);
      }
    return done;
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
namespace FractalSpace
{
  void mini_solve3(const vector<Point*>& found,const double& gc)
  {
    static const double dI204=1.0/204.0;
    Point* pa=0;
    Point* pb=0;
    Point* pc=0;
    int nn=0;
    for(auto p : found)
      {
	int count=0;
	for(int ni=0;ni<6;ni++)
	  {
	    Point* pn=p->get_point_ud_0(ni);
	    for(auto pp : found)
	      if(pn == pp)
		count++;
	  }
	if(count == 2)
	  {
	    pb=p;
	    if(nn=0)
	      {
		pa=found[1];
		pc=found[2];
	      }
	    else if(nn=1)
	      {
		pa=found[0];
		pc=found[2];
	      }
	    else
	      {
		pa=found[0];
		pc=found[1];
	      }
	    break;
	  }
	nn++;
      }
    double suma=-pb->get_potential_point();
    double sumb=-pa->get_potential_point()-pc->get_potential_point();
    double sumc=suma;
    for(int ni=0;ni<6;ni++)
      {
	suma+=pa->get_point_ud_0(ni)->get_potential_point();
	sumb+=pb->get_point_ud_0(ni)->get_potential_point();
	sumc+=pc->get_point_ud_0(ni)->get_potential_point();
      }
    double AB=suma-gc*pa->get_density_point();
    double CD=sumb-gc*pb->get_density_point();
    double EF=sumc-gc*pc->get_density_point();
    pa->set_potential_point((AB*35.0+CD*6.0+EF)*dI204);
    pb->set_potential_point((AB*6.0+CD*36.0+EF*6.0)*dI204);
    pc->set_potential_point((AB+CD*6.0+EF*35.0)*dI204);
  }
}
