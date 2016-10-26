#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void go_ahead_points(vector <Point*>& adj,vector<bool>& ins,vector<bool>& go_ahead)
  {
    //--------------------------------------------------------------------------------------------------------------------------------
    // Given the above information and the provisional list of high_point neighbors, decide which points to generate
    //--------------------------------------------------------------------------------------------------------------------------------
    ins.assign(27,false);
    go_ahead.assign(27,true);
    for(int p_l=0;p_l < 27;++p_l)
      {
	//	ins[p_l]=false;
	//	go_ahead[p_l]=true;
	for(vector <int>::const_iterator ni=Point::nextt[p_l].begin();ni != Point::nextt[p_l].end();ni++)
	  {
	    go_ahead[p_l]= adj[*ni] == 0;
	    if(!go_ahead[p_l]) break;
	  }
      }
    //--------------------------------------------------------------------------------------------------------------------------------
    // A point is inside only if the point is completely covered by neighbor high_cells
    //--------------------------------------------------------------------------------------------------------------------------------
    ins[4]=adj[4];
    ins[10]=adj[10];
    ins[12]=adj[12];
    ins[13]=true;
    ins[1]=adj[1] && adj[4] && adj[10];
    ins[3]=adj[3] && adj[4] && adj[12];
    ins[9]=adj[9] && adj[10] && adj[12];
    ins[0]=adj[0] && ins[1] && ins[3] && ins[9];    
  }
}
