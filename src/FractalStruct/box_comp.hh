#ifndef _Box_Comp_Defined_
#define _Box_Comp_Defined_
namespace FractalSpace
{
  struct box_comp
  {
    bool operator()(const array<int,6>& Boxa,const array<int,6>& Boxb) const
    {
      bool OK=Boxa[2]-Boxb[2] < 0;
      if(!OK)
	{
	  OK=Boxa[1]-Boxb[1] < 0;
	  if(!OK)
	    OK=Boxa[0]-Boxb[0] < 0;
	}
      if(!OK)
	return false;
      return (Boxa[0] > Boxb[1] || Boxa[1] < Boxb[0] ||
	      Boxa[2] > Boxb[3] || Boxa[3] < Boxb[2] ||
	      Boxa[4] > Boxb[5] || Boxa[5] < Boxb[4]);
    }
  };
}
#endif
