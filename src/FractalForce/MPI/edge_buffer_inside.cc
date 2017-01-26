#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void edge_buffer_inside(vector <int>& n,bool& inside,vector <bool>& buff,vector <bool>& edge)
  {
    //MPIrun??
    for(int coord=0;coord<3;coord++)
      {
	buff[coord]=false;
	edge[coord]=false;
	int coord2=coord*2;
	int coord21=coord2+1;
	if(Buffer[coord2]!=0)
	  {
	    if(n[coord] == BBox[coord2])
	      buff[coord2]=true;
	    else if(n[coord] == Box[coord2])
	      edge[coord2]=true;
	  }
	if(Buffer[coord21]!=0)
	  {
	    if(n[coord] == BBox[coord21])
	      buff[coord21]=true;
	    else if(n[coord] == Box[coord21])
	      edge[coord21]=true;
	  }
      }
    for(int coord=0;coord<6;coord++)
      {
	if(edge[coord] && (buff[0] || buff[1] || buff[2]|| buff[3] || buff[4] || buff[5]))
	  edge[coord]=false;
      }
    vector <bool>insides(6,false);
    for(int coord=0;coord<3;coord++)
      {
	if(periods[coord] || (n[coord] > BBox[2*coord] && n[coord] < BBox[2*coord+1]))
	  insides[coord]=true;
      }
    inside=insides[0] && insides[1] && insides[2];
  }
}
