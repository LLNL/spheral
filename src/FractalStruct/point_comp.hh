#ifndef _Point_Comp_Defined_
#define _Point_Comp_Defined_
namespace FractalSpace
{
  typedef std::pair<std::array<int,3>,Point*> pap;
  struct point_comp2
  {
    bool operator()(const std::array<int,3>& pa,const std::array<int,3>& pb) const
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
    bool operator()(const std::array<int,3>& pa,const std::array<int,3>& pb) const
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
    bool operator()(const std::array<int,3>& pa,const std::array<int,3>& pb) const
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
  struct point_comp4
  {
    bool operator()(const std::array<int,4>& pa,const std::array<int,4>& pb) const
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
  template <int N>  struct point_compN
  {
    bool operator()(const std::array<int,N>& pa,const std::array<int,N>& pb) const
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
  struct vector_comp_up
  {
    template <class T> bool operator()(const std::vector<T*> vA,const std::vector<T*> vB) const
    {
      return vA.size() < vB.size();
    }
  };
  struct vector_comp_down
  {
    template <class T> bool operator()(const std::vector<T*> vA,const std::vector<T*> vB) const
    {
      return -vA.size() < -vB.size();
    }
  };
}
#endif
