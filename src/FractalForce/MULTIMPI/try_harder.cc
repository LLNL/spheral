#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  Point* try_harder(Point& point0,const int& ni, const bool& easy)
  {
    Point* adj1=0;
    int rp0=point0.get_real_pointer();
    int spam=rp0*27+ni;
    int spam8=spam*8;
    int tries=Point::dupes[spam];
    Point* pl0=&point0;
    if(!easy && !Point::left[rp0]) 
      pl0=pl0->move_rp(13);
    pl0=pl0->move_rp(0);
      //
    Point* ph0=pl0->get_point_pointer();
    Point* ph=0;
    vector <bool> eureka_h(27);
    ph0->get_eureka_adj(eureka_h);
    for(int t=0;t<tries;t++)
      {
	int witch=Point::phl[spam8+t];
	int rh=witch/27;
	if(!eureka_h[rh])
	  continue;
	ph=ph0->move_adj(0,Point::corner_a[rh]);
	ph=ph->move_adj(Point::corner_b[rh],0);
	Point* pl=ph->get_p_daughter_point();
	if(pl->get_real_pointer() != 0) continue;
	int rl=witch % 27;
	if(!ph->get_eureka_dau(rl)) continue;
	adj1=pl;
	if(!easy && !Point::left[rl])
	  adj1=adj1->move_rp(13);
	adj1=adj1->move_rp(rl);
	return adj1;
      }
    return 0;
  }
}
