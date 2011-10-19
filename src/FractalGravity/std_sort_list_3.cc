namespace std_sort
{
  int dim0=0;
  int dim1=1;
  int dim2=2;
}
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void std_sort_list_3(Fractal& fractal,Group& group,const bool& remove_duplicates,const int& what)
  {
    using namespace std_sort;
    vector <int>p(3,-1);
    vector <int>pn(3,-1);
    bool period=fractal.get_periodic();
    int skip=Misc::pow(2,fractal.get_level_max()-group.get_level());
    int total=fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max());
    //  cout << "total= " << total << endl;
    if(what == 0)
      {
	std_sort::dim0=0;
	std_sort::dim1=1;
	std_sort::dim2=2;
      }
    else if(what == 1)
      {
	std_sort::dim0=1;
	std_sort::dim1=2;
	std_sort::dim2=0;
      }
    else if(what == 2)
      {
	std_sort::dim0=2;
	std_sort::dim1=0;
	std_sort::dim2=1;
      }
    else
      {
	assert(0);
      }
    int dim0t2p1=dim0*2+1;
    group.list_points.sort(LesserPoint);
    if(remove_duplicates)
      {
	list <Point*>::iterator point_itr = group.list_points.begin();
	while(point_itr != group.list_points.end())
	  {
	    if((*point_itr)->get_dupe())
	      point_itr=group.list_points.erase(point_itr);
	    else
	      point_itr++;
	  }
      }
    list <Point*>::iterator final=group.list_points.end();
    final--;
    Point* left=0;
    for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr!=final;point_itr++)
      {
	list <Point*>::const_iterator next_itr=point_itr;
	next_itr++;
	(*point_itr)->get_pos_point(p);
	(*next_itr)->get_pos_point(pn);
	//      cout << "std a " << p[0] << " " << p[1] << " " << p[2] << " " << pn[0] << " " << pn[1] << " " << pn[2] << " " << skip << " " << (*point_itr)->get_real_pointer() << " " << *point_itr <<  endl;
	if(p[dim0] == 0) left=*point_itr;
	//      if(left) cout << "std a left " << px << " " << py << " " << pz << endl;
	if(p[dim2] == pn[dim2] && p[dim1] == pn[dim1])
	  {
	    if(p[dim0]+skip == pn[dim0])
	      (*point_itr)->set_point_ud(*next_itr,dim0t2p1);
	  }
	if(period && pn[dim0]+skip == total)
	  {
	    if(left && left->get_pos_point(dim1) == pn[dim1] && left->get_pos_point(dim2) == pn[dim2])
	      (*next_itr)->set_point_ud(left,dim0t2p1);
	  }
      }
  }
  bool LesserPoint(Point* p1,Point* p2)
  {
    using namespace std_sort;
    int dz=p1->get_pos_point(std_sort::dim2)-p2->get_pos_point(std_sort::dim2);
    if(dz < 0 )
      return true;
    else if(dz > 0)
      return false;
    else
      {
	int dy=p1->get_pos_point(std_sort::dim1)-p2->get_pos_point(std_sort::dim1);
	if(dy < 0 )
	  return true;
	else if(dy > 0)
	  return false;
	else
	  {
	    int dx=p1->get_pos_point(std_sort::dim0)-p2->get_pos_point(std_sort::dim0);
	    if(dx < 0 )
	      return true;
	    else if(dx > 0)
	      return false;
	    else
	      if(p1->get_real_pointer() < p2->get_real_pointer())
		{
		  p2->set_dupe(true);
		  return true;
		}
	      else if(p1->get_real_pointer() > p2->get_real_pointer())
		{
		  p1->set_dupe(true);
		  return false;
		}
	      else
		{
		  cout << "bad std_sort " << p1 << " " << p1->get_real_pointer() << " " << p2 << " " << p2->get_real_pointer() ;
		  cout << " " << dx << " " << dy << " " << dz << endl;
		  assert(0);
		}
	  }
      }   
    return true;
  }
}
