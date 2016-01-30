#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_clever_boxes(Group& group,vector <Point*>p_points_left,vector <Point*>p_points_right)
  {
    counst int lev=group.get_level();
    vector <Point*> listP;
    for(vector <Point*>::iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* P=(*point_itr)->get_point_pointer();
	if(P != 0)
	  listP.push_back(P);
      }
    sort3_list(listP,0);

    vector <int> real_left(7);
    vector <int> real_right(7);
    vector <bool> belongs_to_me(27);
    vector <Point*> pointers(27);
    for(vector <Point*>::iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	if(p_point->get_real_pointer() != 0)
	  continue;
	bool xwall=false;
	bool ywall=false;
	p_point->all_mine(pointers,belongs_to_me);
	int boxes=0;
	real_left[boxes]=0;
	real_right[boxes]=13;
	boxes++;
	if(belongs_to_me[2])
	  {
	    real_left[boxes]=2;
	    real_right[boxes]=14;
	    xwall=true;
	    if(belongs_to_me[8])
	      real_right[boxes]=17;
	    boxes++;
	  }
	if(belongs_to_me[6])
	  {
	    real_left[boxes]=6;
	    real_right[boxes]=16;
	    ywall=true;
	    if(!xwall && belongs_to_me[8])
	      real_right[boxes]=17;
	    boxes++;
	  }
	if(!xwall && !ywall && belongs_to_me[8])
	  {
	    real_left[boxes]=8;
	    real_right[boxes]=17;
	    boxes++;
	  }
	xwall=false;
	ywall=false;
	if(belongs_to_me[18])
	  {
	    real_left[boxes]=18;
	    real_right[boxes]=22;
	    boxes++;
	  }
	if(belongs_to_me[20])
	  {
	    real_left[boxes]=20;
	    real_right[boxes]=23;
	    xwall=true;
	    if(belongs_to_me[26])
	      real_right[boxes]=26;
	    boxes++;
	  }
	if(belongs_to_me[24])
	  {
	    real_left[boxes]=24;
	    real_right[boxes]=25;
	    ywall=true;
	    if(!xwall && belongs_to_me[26])
	      real_right[boxes]=26;
	    boxes++;
	  }
	if(!xwall && !ywall && belongs_to_me[26])
	  {
	    //	    cout << "this cannot be right" << endl;
	    assert(0);
	  }
	for(int ni=0;ni<boxes;ni++)
	  {
	    p_points_left.push_back(pointers[real_left[ni]]);
	    p_points_right.push_back(pointers[real_right[ni]]);
	  }
      }
  }
}
