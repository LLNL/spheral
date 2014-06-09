#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool test_up_down(Group& group)
  {
    //  vector <Point*> ud(6);
    bool badx=false;
    bool bady=false;
    bool badz=false;
    for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	int rp=p_point->get_real_pointer();
	//      p_point->get_point_ud(ud);
	int rpx=rp % 3;
	int rpy=rp/3 % 3;
	int rpz=rp/9;
	badx= rpx < 2 && p_point->get_point_ud(1) == 0;
	bady= rpy < 2 && p_point->get_point_ud(3) == 0;
	badz= rpz < 2 && p_point->get_point_ud(5) == 0;
	if(badx || bady || badz) 
	  {
	    cout << "baaad uppers " << endl;
	    cout << p_point << endl;
	    cout << badx << bady << badz << endl;
	    cout << rp << " " << rpx << " " << rpy << " " << rpz << endl;
	    cout << p_point->get_point_ud(1) << endl;
	    cout << p_point->get_point_ud(3) << endl;
	    cout << p_point->get_point_ud(5) << endl; 
	    p_point->dump();
	    return true;
	  }
	badx= rpx > 0 && p_point->get_point_ud(0) == 0;
	bady= rpy > 0 && p_point->get_point_ud(2) == 0;
	badz= rpz > 0 && p_point->get_point_ud(4) == 0;
	if(badx || bady || badz) 
	  {
	    cout << "baaad downers " << endl;
	    cout << p_point << endl;
	    cout << badx << bady << badz << endl;
	    cout << rp << " " << rpx << " " << rpy << " " << rpz << endl;
	    cout << p_point->get_point_ud(0) << endl;
	    cout << p_point->get_point_ud(2) << endl;
	    cout << p_point->get_point_ud(4) << endl; 
	    p_point->dump();
	    return true;
	  }
	badx=rp == 14 && p_point->get_point_ud(1) != 0;
	bady=rp == 16 && p_point->get_point_ud(3) != 0;
	badz=rp == 22 && p_point->get_point_ud(5) != 0;
	if(badx || bady || badz) 
	  {
	    cout << "baaad faces " << endl;
	    cout << p_point << endl;
	    cout << badx << bady << badz << endl;
	    p_point->dump();
	    return true;
	  }
	//      cout << badx << bady << badz << endl;
	//      p_point->dump();
      }
    return false;
  }
}
