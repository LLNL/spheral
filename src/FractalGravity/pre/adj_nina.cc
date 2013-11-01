#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void adj_nina(Point& point,vector <Point*>& adj)
  {
    adj[0]=point.neighbor(0,2,4);
    adj[1]=point.neighbor(2,4);
    adj[2]=point.neighbor(1,2,4);
    adj[3]=point.neighbor(0,4);
    adj[4]=point.neighbor(4);
    adj[5]=point.neighbor(1,4);
    adj[6]=point.neighbor(0,3,4);
    adj[7]=point.neighbor(3,4);
    adj[8]=point.neighbor(1,3,4);
    adj[9]=point.neighbor(0,2);
    adj[10]=point.neighbor(2);
    adj[11]=point.neighbor(1,2);
    adj[12]=point.neighbor(0);
    adj[13]=&point;
    adj[14]=point.neighbor(1);
    adj[15]=point.neighbor(0,3);
    adj[16]=point.neighbor(3);
    adj[17]=point.neighbor(1,3);
    adj[18]=point.neighbor(0,2,5);
    adj[19]=point.neighbor(2,5);
    adj[20]=point.neighbor(1,2,5);
    adj[21]=point.neighbor(0,5);
    adj[22]=point.neighbor(5);
    adj[23]=point.neighbor(1,5);
    adj[24]=point.neighbor(0,3,5);
    adj[25]=point.neighbor(3,5);
    adj[26]=point.neighbor(1,3,5);
    //
    vector <bool> sure(27);
    for(int ni=0;ni<27;ni++)
      sure[ni]=adj[ni] != 0;
    //
    sure[4]=true;
    sure[22]=true;
    sure[10]=true;
    sure[16]=true;
    sure[12]=true;
    sure[14]=true;
    //
    sure[9]= sure[9] || adj[10] != 0 || adj[12] != 0;
    sure[11]= sure[11] || adj[10] != 0 || adj[14] != 0;
    sure[15]= sure[15] || adj[12] != 0 || adj[16] != 0;
    sure[17]= sure[17] || adj[14] != 0 || adj[16] != 0;
    sure[3]= sure[3] || adj[4] != 0 || adj[12] != 0;
    sure[5]= sure[5] || adj[4] != 0 || adj[14] != 0;
    sure[21]= sure[21] || adj[12] != 0 || adj[22] != 0;
    sure[23]= sure[23] || adj[14] != 0 || adj[22] != 0;
    sure[1]= sure[1] || adj[4] != 0 || adj[10] != 0;
    sure[7]= sure[7] || adj[4] != 0 || adj[16] != 0;
    sure[19]= sure[19] || adj[10] != 0 || adj[22] != 0;
    sure[25]= sure[25] || adj[16] != 0 || adj[22] != 0;
    //
    sure[0]= sure[0] || adj[1] != 0 || adj[3] != 0 || adj[9] != 0;
    sure[2]= sure[2] || adj[1] != 0 || adj[5] != 0 || adj[11] != 0;
    sure[6]= sure[6] || adj[7] != 0 || adj[3] != 0 || adj[15] != 0;
    sure[8]= sure[8] || adj[7] != 0 || adj[5] != 0 || adj[17] != 0;
    sure[18]= sure[18] || adj[19] != 0 || adj[21] != 0 || adj[9] != 0;
    sure[20]= sure[20] || adj[19] != 0 || adj[23] != 0 || adj[11] != 0;
    sure[24]= sure[24] || adj[25] != 0 || adj[21] != 0 || adj[15] != 0;
    sure[26]= sure[26] || adj[25] != 0 || adj[23] != 0 || adj[17] != 0;
    //
    Group* h_group0=point.get_p_in_high_group();
    vector <bool>eureka(27,false);
    for(int ni=0;ni<27;ni++)
      {
	if(adj[ni] != 0 && !(adj[ni]->get_it_is_high() && h_group0 == adj[ni]->get_p_in_high_group()))
	  adj[ni]=0;
	if(!sure[ni])
	  {
	    adj[ni]=try_harder(point,ni,true);
	    if(adj[ni] != 0 && !(adj[ni]->get_it_is_high() && h_group0 == adj[ni]->get_p_in_high_group()))
	      adj[ni]=0;
	    sure[ni]=true;  
	  }
	  eureka[ni]=adj[ni] != 0;
      }
    point.set_eureka_adj(eureka);
  }
}
