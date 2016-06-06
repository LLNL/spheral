#ifndef _Point_Comp_Defined_
#define _Point_Comp_Defined_
namespace FractalSpace
{
  typedef std::pair<array<int,3>,Point*> pap;
  struct point_comp2
  {
    bool operator()(const array<int,3>& pa,const array<int,3>& pb) const
    {
      int dif=pa[2]-pb[2];
      if(dif != 0)
	return dif < 0;
      dif=pa[1]-pb[1];
      if(dif != 0)
	return dif < 0;
      return pa[0]-pb[0] < 0;
    }
  };
  struct point_comp1
  {
    bool operator()(const array<int,3>& pa,const array<int,3>& pb) const
    {
      int dif=pa[1]-pb[1];
      if(dif != 0)
	return dif < 0;
      dif=pa[0]-pb[0];
      if(dif != 0)
	return dif < 0;
      return pa[2]-pb[2] < 0;
    }
  };
  struct point_comp0
  {
    bool operator()(const array<int,3>& pa,const array<int,3>& pb) const
    {
      int dif=pa[0]-pb[0];
      if(dif != 0)
	return dif < 0;
      dif=pa[2]-pb[2];
      if(dif != 0)
	return dif < 0;
      return pa[1]-pb[1] < 0;
    }
  };
  struct pointer_comp
  {
    bool operator()(const Point* Pa,const Point* Pb) const
    {
      int dif=Pa->get_pos_point_z()-Pb->get_pos_point_z();
      if(dif != 0)
	return dif < 0;
      dif=Pa->get_pos_point_y()-Pb->get_pos_point_y();
      if(dif != 0)
	return dif < 0;
      return Pa->get_pos_point_x()-Pb->get_pos_point_x() < 0;
    }
  };
}
#endif
