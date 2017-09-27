#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void assign_buffer_edge(Group& group,Fractal& fractal)
  {
    int lev=group.get_level();
    int lmax=fractal.get_level_max();
    int zoom=Misc::pow(2,lmax);
    int delta=Misc::pow(2,lmax-lev);
    vector <int> Buffer(6);
    fractal.getBuffer(Buffer);
    vector <int> Box(6);
    fractal.getBox(Box);
    vector <int> BBox(6);
    fractal.getBBox(BBox);
    vector <int>pos(3);
    vector <int>lefte(3);
    vector <int>righte(3);
    vector <int>leftb(3);
    vector <int>rightb(3);
    vector <bool>testl(3);
    vector <bool>testr(3);
    bool pass=false;
    vector <bool>buff(6,false);
    vector <bool>edge(6,false);
    for(int n=0;n<3;n++)
      {
	testl[n]=Buffer[n*2] > 0;
	testr[n]=Buffer[n*2+1] > 0;
	lefte[n]=Box[n*2]*zoom;
	leftb[n]=lefte[n]-delta;
	rightb[n]=BBox[n*2+1]*zoom;
	righte[n]=rightb[n]-delta;
      }
    bool bg=false;
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	point.get_pos_point(pos);
	pass=false;
	buff.assign(6,false);
	edge.assign(6,false);
	for(int n=0;n<3;n++)
	  {
	    int n2=n*2;
	    if(testl[n])
	      {
		pass=pass || pos[n] < leftb[n];
		buff[n2]=pos[n] == leftb[n];
		edge[n2]=pos[n] == lefte[n];
		bg=bg || pos[n] <= lefte[n];
	      }
	    if(testr[n])
	      {
		pass=pass || pos[n] > rightb[n];
		buff[n2+1]=pos[n] == rightb[n];
		edge[n2+1]=pos[n] == righte[n];
		bg=bg || pos[n] >= righte[n];
	      }
	  }
	for(int n=0;n<6;n++)
	  {
	    edge[n]=edge[n] && !(buff[0] || buff[1] || buff[2] || buff[3] || buff[4] || buff[5] || pass);
	    buff[n]=buff[n] && !pass;
	  }
	point.set_passive_point(pass);
	point.set_buffer_point(buff);
	point.set_edge_point(edge);
      }    
    group.set_buffer_group(bg);
  }
}
