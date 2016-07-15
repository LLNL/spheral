#ifndef _Box_Defined_
#define _Box_Defined_
namespace FractalSpace
{
  class Box
  {
    vector <double> boxC;
    vector <double> boxW;
  public:
    list <Box*>overlaps;
    Box()
    {
      boxC.assign(3,0.5);
      boxW.assign(3,1.0);
    }
    Box(const vector <double>& C,const vector <double>& W)
    {
      boxC=C;
      boxW=W;
    }
     void setC(const vector <double>& C)
    {
      boxC=C;
    }
     void setW(const vector <double>& W)
    {
      boxW=W;
    }
     void getC(vector <double>& C)
    {
      C=boxC;
    }
     void getW(vector <double>& W)
    {
      W=boxW;
    }
    bool overlap(Box& boxx, const int& dir)
    {
      return 0.5*(boxW[dir]+boxx.boxW[dir]) >= abs(boxC[dir]-boxx.boxC[dir]);
    }
    bool overlap(Box& boxx)
    {
      if(0.5*(boxW[0]+boxx.boxW[0]) < abs(boxC[0]-boxx.boxC[0]))
	{
	  if(0.5*(boxW[1]+boxx.boxW[1]) < abs(boxC[1]-boxx.boxC[1]))
	    {
	      if(0.5*(boxW[2]+boxx.boxW[2]) < abs(boxC[2]-boxx.boxC[2]))
		{
		  return false;	  
		}
	    }
	}
      return true;
    }
    ~Box()
    {}
  };
}
#endif
/*
#ifndef _Cube_Defined_
#define _Cube_Defined_
namespace FractalSpace
{
  class Cube : protected Box
  {
    vector <double> cubeC;
    double cubeW;
  public:
    Cube()
    {
      cubeC.assign(3,0.5);
      cubeW=1.0;
    }
    Cube(const vector <double>& C,const double& W)
    {
      cubeC=C;
      cubeW=W;
    }
    ~Cube()
    {}
  };
}
#endif
*/
