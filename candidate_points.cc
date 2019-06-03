#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool Point::calc_candidates=true;
  vector <bool> Point::left(27,false);
  vector <bool> Point::corner(27,false);
  vector <bool> Point::edge(27,false);
  vector <bool> Point::face(27,false);
  vector <bool> Point::center(27,false);
  vector < vector <int> > Point::cefc(27);
  vector <int>Point::phl(5832,0);
  vector <int>Point::dupes(729,0);
  vector <int>Point::corner_a;
  vector <int>Point::corner_b;
  vector <int>Point::sequence;
  vector <int>Point::updown;
  vector<vector<int>> Point::order;
  vector <vector <int> >Point::nextt(27);

  void candidate_points()
  {
    for(int niz=0;niz<2;niz++)
      for(int niy=0;niy<2;niy++)
	for(int nix=0;nix<2;nix++)
	  Point::left[nix+niy*3+niz*9]=true;

    Point::cefc[1]={0,1};
    Point::cefc[3]={2,3};
    Point::cefc[4]={0,1,2,3};
    Point::cefc[5]={2,3};
    Point::cefc[7]={0,1};
    Point::cefc[9]={4,5};
    Point::cefc[10]={0,1,4,5};
    Point::cefc[11]={4,5};
    Point::cefc[12]={2,3,4,5};
    Point::cefc[14]={2,3,4,5};
    Point::cefc[15]={4,5};
    Point::cefc[16]={0,1,4,5};
    Point::cefc[17]={4,5};
    Point::cefc[19]={0,1};
    Point::cefc[21]={2,3};
    Point::cefc[22]={0,1,2,3};
    Point::cefc[23]={2,3};
    Point::cefc[25]={0,1};
    
    Point::corner[0]=true;
    Point::edge[1]=true;
    Point::corner[2]=true;
    Point::edge[3]=true;
    Point::face[4]=true;
    Point::edge[5]=true;
    Point::corner[6]=true;
    Point::edge[7]=true;
    Point::corner[8]=true;
    Point::edge[9]=true;
    Point::face[10]=true;
    Point::edge[11]=true;
    Point::face[12]=true;
    Point::center[13]=true;
    Point::face[14]=true;
    Point::edge[15]=true;
    Point::face[16]=true;
    Point::edge[17]=true;
    Point::corner[18]=true;
    Point::edge[19]=true;
    Point::corner[20]=true;
    Point::edge[21]=true;
    Point::face[22]=true;
    Point::edge[23]=true;
    Point::corner[24]=true;
    Point::edge[25]=true;
    Point::corner[26]=true;
    
    Point::order.push_back(vector<int>());
    Point::order.back()={0,2,1,4,0,3,1};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={1,2,0,4,1,3,0};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={0,3,1,4,0,2,1};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={1,3,0,4,1,2,0};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={0,2,1,5,0,3,1};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={1,2,0,5,1,3,0};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={0,3,1,5,0,2,1};
    
    Point::order.push_back(vector<int>());
    Point::order.back()={1,3,0,5,1,2,0};
    
    Point::corner_a={0,0,1,0,0,1,2,2,3,0,0,1,0,0,1,2,2,3,4,4,5,4,4,5,6,6,7};
    
    Point::corner_b={7,6,6,5,4,4,5,4,4,3,2,2,1,0,0,1,0,0,3,2,2,1,0,0,1,0,0};
    
    Point::sequence={0,1,3,2,6,7,5,4};
    Point::updown={1,3,0,5,1,2,0,4};
    
    for(int pl=0;pl < 27;pl++)
      {
	int plz=pl/9;
	int ply=pl/3 %3;
	int plx=pl % 3;
	for(int pla=0;pla < 27;pla++)
	  {
	    int plaz=pla/9-1;
	    int play=pla/3 %3 -1;
	    int plax=pla % 3 -1;
	    int z=plz+plaz;
	    int y=ply+play;
	    int x=plx+plax;
	    int xyz=x+3*y+9*z;
	    for(int ph=0;ph < 27;ph++)
	      {
		int phz=ph/9;
		int phy=ph/3 %3;
		int phx=ph % 3;
		for(int plb=0;plb < 27;plb++)
		  {
		    int plbz=plb/9;
		    int plby=plb/3 %3;
		    int plbx=plb % 3;
		    int zb=plbz+(phz-1)*2;
		    int yb=plby+(phy-1)*2;
		    int xb=plbx+(phx-1)*2;
		    if(z == zb && y==yb && x==xb)
		      {
			bool ins=abs(z-1) <=1 && abs(y-1) <=1 && abs(x-1) <=1;
			if(! ins ||(ins && plb <= xyz))
			  {
			    int which=(pl*27+pla)*8;
			    Point::phl[which+Point::dupes[pl*27+pla]]=ph*27+plb;
			    Point::dupes[pl*27+pla]++;
			  }
		      }
		  }
	      }
	  }
      }
    array<array<int,27>,27> positions;
    array<array<bool,27>,27> decisions;
    //--------------------------------------------------------------------------------------------------------------------------------
    // Find which of the 27 points in a 3x3 cube can be affected by the 3x3 cube of high points.
    // need only do this once.
    //--------------------------------------------------------------------------------------------------------------------------------
    for(int p_z_h=0;p_z_h < 3;++p_z_h)
      {
	int pos_z_h=(p_z_h-1)*2;
	for(int p_y_h=0;p_y_h < 3;++p_y_h)
	  {
	    int pos_y_h=(p_y_h-1)*2;
	    for(int p_x_h=0;p_x_h < 3;++p_x_h)
	      {
		int pos_x_h=(p_x_h-1)*2;
		int pos_h=p_z_h*9+p_y_h*3+p_x_h;
		for(int p_z_l=0;p_z_l < 3;++p_z_l)
		  {
		    int pos_z_l=p_z_l;
		    int d_z=pos_z_l-pos_z_h;
		    for(int p_y_l=0;p_y_l < 3;++p_y_l)
		      {
			int pos_y_l=p_y_l;
			int d_y=pos_y_l-pos_y_h;
			for(int p_x_l=0;p_x_l < 3;++p_x_l)
			  {
			    int pos_x_l=p_x_l;
			    int d_x=pos_x_l-pos_x_h;
			    int pos_l=p_z_l*9+p_y_l*3+p_x_l;
			    positions[pos_l][pos_h]=-1;
			    decisions[pos_l][pos_h]=
			      abs(d_x-1) <=1 && 
			      abs(d_y-1) <=1 && 
			      abs(d_z-1) <=1 ;
			    if(decisions[pos_l][pos_h])
			      positions[pos_l][pos_h]=d_x+d_y*3+d_z*9;
			  }
		      }
		  }
	      }
	  }
      }
    for (int p_l=0;p_l < 27;++p_l)
      {
	for(int p_h=0;p_h < 27; ++p_h)
	  {
	    if(p_h == 13) continue;
	    if(decisions[p_l][p_h]  && positions[p_l][p_h] < p_l)
	      Point::nextt[p_l].push_back(p_h);
	  }
      }
    Point::calc_candidates=false;
  }
}
