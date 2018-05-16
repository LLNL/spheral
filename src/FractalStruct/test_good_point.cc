#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool test_good_point(Point* p1,Fractal_Memory& mem,int level)
  {
    bool good=p1->get_inside();
    if(good)
      {
	vector <int>pos(3);
	p1->get_pos_point(pos);
	good=vector_in_box(pos,mem.BoxesLev[mem.p_mess->FractalRank][level]);
	if(!good)
	  {
	    for(auto FR : mem.Touchy)
	      {
		good=vector_in_box(pos,mem.BoxesLev[FR][level]);
		if(good)
		  break;
	      }
	    assert(good);
	  }
      }
    return good;
  }
}
