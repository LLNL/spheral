#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_big_boxes(Fractal& frac,Group& group,vector <Point*>p_points_left,vector <Point*>p_points_right)
  {
    const int level=group.get_level();
    const int spacing=Misc::pow(2,fractal.get_level_max()-level);
    vector <Point*> pleft0;
    vector <Point*> pright0;
    vector <Point*> pleft1;
    vector <Point*> pright1;
    for(vector <Point*>::iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p->get_real_pointer() != 12 || p->get_point_ud(0) != 0)
	  continue;
	pleft0.push_back(p->get_point_ud_0[2])->get_point_ud_0[4];
      }
    int psize0=pleft0.size();
    pright0.resize(psize0);
    for(int ni=0;ni < psize0;ni++)
      {
	Point* p=(pleft0[ni]->get_point_ud_0[3])->get_point_ud_0[5];
	while(p->get_real_pointer() != 14)
	  p=(p->get_point_ud_0[1])->get_point_ud_0[1];
	pright0[ni]=p;
      }
    vector <bool>edge_y;
    for(int ni=0;ni < psize0;ni++)
      {
	edge_y.clear();
	Point* p=pleft0[ni]->move_rp(16);
	bool more=true;
	while(more)
	  {
	    edge_y.push_back(p->get_point_ud(3)==0)
	    p=(p->get_point_ud_0(1))->get_point_ud(1);
	    more=p != 0;
	  }
	int psize1=edge_y.size();

	Point* p0=pleft0[ni];
	pleft1.push_back(p1);
	int rp=13;
	if(edge_y[0])
	  rp=16;
	Point* p1=p0.move_rp(rp);
	
	

      }
  }
}
