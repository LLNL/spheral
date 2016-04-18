#ifndef _Point_Comp_Defined_
#define _Point_Comp_Defined_
namespace FractalSpace
{
  typedef std::pair<array<int,3>,Point*> pap;
  struct point_comp
  {
    bool operator()(const array<int,3>& pa,const array<int,3>& pb) const
    {
      int dif=pa[2]-pb[2];
      if(dif != 0)
	return dif < 0;
      dif=pa[1]-pb[1];
      if(dif != 0)
	return dif < 0;
      return pa[1]-pb[0] < 0;
    }
  };
}
#endif
