#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool test_group(Group& group)
  {
    bool badd=false;
    vector <Point*> ud(6);
    // for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
    for(auto pp : group.list_points)
      {
	// Point& point=**point_itr;
	const int rp=pp->get_real_pointer();
	pp->get_point_ud(ud);
	int sumn=0;
	int ni0=0;
	int ni1=1;
	bool badd1=false;
	for(int ni=0;ni<3;ni++)
	  {
	    if(ud[ni0] != 0) sumn++;
	    if(ud[ni1] != 0) sumn++;
	    badd1= badd1 ||
	      (ud[ni0] != 0 && ud[ni0]->get_point_ud(ni1) != pp ) ||
	      (ud[ni1] != 0 && ud[ni1]->get_point_ud(ni0) != pp );
	    ni0+=2;
	    ni1+=2;
	  }
	int rpx=rp % 3;
	int rpy=(rp/3) % 3;
	int rpz=rp/9;
	int sumr=(rpx % 2) + (rpy % 2) + (rpz % 2);
	bool badd2=(sumr == 0 && sumn < 3) || (sumr == 1 && sumn < 4) || 
	  (sumr == 2 && sumn < 5) || (sumr == 3 && sumn < 6);
	bool badd3=(rp==14 && ud[1] != 0) || (rp==16 && ud[3] != 0) || (rp==22 && ud[5] != 0);
	bool badd4=pp->get_inside() && sumn < 6;
	bool badd5= (sumr == 0 && pp->get_point_pointer_t() == 0) ||(sumr != 0 && pp->get_point_pointer_t() != 0);
	bool badd6=(rpx < 2 && ud[1] == 0) ||(rpy < 2 && ud[3] == 0) ||(rpz < 2 && ud[5] == 0) ||
	  (rpx > 0 && ud[0] == 0) ||(rpy > 0 && ud[2] == 0) ||(rpz > 0 && ud[4] == 0);
	badd=badd1 || badd2 || badd3 || badd4 || badd5 || badd6;
	if(badd)
	  {
	    //	    cerr << "badd " << badd1 << badd2 << badd3 << badd4 << badd5 << badd6 << "\n";
	    pp->dump();
	    assert(!badd);
	  }
      }
    return badd;
  }
}
