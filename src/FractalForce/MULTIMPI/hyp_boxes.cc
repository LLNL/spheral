#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_big_boxes(vector <Point*>p_points_left,vector <Point*>p_points_right)
  {
    vector <Point*>tmp_left_small;
    vector <Point*>tmp_right_small;
    vector <Point*>tmp_left_big;
    vector <Point*>tmp_right_big;
    vector <Point*>tmp_left;
    vector <Point*>tmp_right;
    int boxes=p_points_left.size();
    vector <int>order;
    int boxes_big=0;
    for(int box=0;box<boxes;box++)
      {
	if(p_points_left[box]->get_real_pointer()==0)
	  {
	    tmp_left_big.push_back(p_points_left[box]);
	    tmp_right_big.push_back(p_points_right[box]);
	    order.push_back(boxes_big);
	    boxes_big++;
	  }
	else
	  {
	    tmp_left_small.push_back(p_points_left[box]);
	    tmp_right_small.push_back(p_points_right[box]);
	  }
      }
    int dir=0;
    for(int box=0;box<boxes_big-1;box+=2)
      {
	Point* pl0=tmp_left_big[box];
	Point* pr0=tmp_right_big[box];
	Point* pl1=tmp_left_big[box+1];
	Point* pr1=tmp_right_big[box+1];
	ok=check_touch(pl0,pr0,pl1,pr2,dir);
      }
  }
  bool check_touch(Point* pl0,Point* pr0,Point* pl1,Point* pr1,const int& dir)
  {
    
  }
}
