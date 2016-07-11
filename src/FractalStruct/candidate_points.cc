#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  //  ofstream& FileFractal=fractal.p_file->FileFractal;
  vector <bool> Point::left(27);
  vector <bool> Point::corner(27);
  vector <bool> Point::edge(27);
  vector <bool> Point::face(27);
  vector <bool> Point::center(27);
  vector < vector <int> > Point::cefc(27);
  vector <int>Point::phl;
  vector <int>Point::dupes;
  vector <int>Point::corner_a;
  vector <int>Point::corner_b;
  vector <int>Point::sequence;
  vector <int>Point::updown;
  int Point::order [8][7];
  vector <vector <int> >Point::nextt(27);
  ofstream* Point::p_FILE;
  void candidate_points()
  {
    Point::left.assign(27,false);
    for(int niz=0;niz<2;niz++)
      {
	for(int niy=0;niy<2;niy++)
	  {
	    for(int nix=0;nix<2;nix++)
	      Point::left[nix+niy*3+niz*9]=true;
	  }
      }
    Point::cefc.resize(27);
    Point::cefc[1].push_back(0);
    Point::cefc[1].push_back(1);
    Point::cefc[3].push_back(2);
    Point::cefc[3].push_back(3);
    Point::cefc[4].push_back(0);
    Point::cefc[4].push_back(1);
    Point::cefc[4].push_back(2);
    Point::cefc[4].push_back(3);
    Point::cefc[5].push_back(2);
    Point::cefc[5].push_back(3);
    Point::cefc[7].push_back(0);
    Point::cefc[7].push_back(1);
    Point::cefc[9].push_back(4);
    Point::cefc[9].push_back(5);
    Point::cefc[10].push_back(0);
    Point::cefc[10].push_back(1);
    Point::cefc[10].push_back(4);
    Point::cefc[10].push_back(5);
    Point::cefc[11].push_back(4);
    Point::cefc[11].push_back(5);
    Point::cefc[12].push_back(2);
    Point::cefc[12].push_back(3);
    Point::cefc[12].push_back(4);
    Point::cefc[12].push_back(5);
    Point::cefc[14].push_back(2);
    Point::cefc[14].push_back(3);
    Point::cefc[14].push_back(4);
    Point::cefc[14].push_back(5);
    Point::cefc[15].push_back(4);
    Point::cefc[15].push_back(5);
    Point::cefc[16].push_back(0);
    Point::cefc[16].push_back(1);
    Point::cefc[16].push_back(4);
    Point::cefc[16].push_back(5);
    Point::cefc[17].push_back(4);
    Point::cefc[17].push_back(5);
    Point::cefc[19].push_back(0);
    Point::cefc[19].push_back(1);
    Point::cefc[21].push_back(2);
    Point::cefc[21].push_back(3);
    Point::cefc[22].push_back(0);
    Point::cefc[22].push_back(1);
    Point::cefc[22].push_back(2);
    Point::cefc[22].push_back(3);
    Point::cefc[23].push_back(2);
    Point::cefc[23].push_back(3);
    Point::cefc[25].push_back(0);
    Point::cefc[25].push_back(1);
    //
    Point::corner.assign(27,false);
    Point::edge.assign(27,false);
    Point::face.assign(27,false);
    Point::center.assign(27,false);
    //
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
    //
    Point::order[0][0]=0;
    Point::order[0][1]=2;
    Point::order[0][2]=1;
    Point::order[0][3]=4;
    Point::order[0][4]=0;
    Point::order[0][5]=3;
    Point::order[0][6]=1;
    //
    Point::order[1][0]=1;
    Point::order[1][1]=2;
    Point::order[1][2]=0;
    Point::order[1][3]=4;
    Point::order[1][4]=1;
    Point::order[1][5]=3;
    Point::order[1][6]=0;
    //
    Point::order[2][0]=0;
    Point::order[2][1]=3;
    Point::order[2][2]=1;
    Point::order[2][3]=4;
    Point::order[2][4]=0;
    Point::order[2][5]=2;
    Point::order[2][6]=1;
    //
    Point::order[3][0]=1;
    Point::order[3][1]=3;
    Point::order[3][2]=0;
    Point::order[3][3]=4;
    Point::order[3][4]=1;
    Point::order[3][5]=2;
    Point::order[3][6]=0;
    //
    Point::order[4][0]=0;
    Point::order[4][1]=2;
    Point::order[4][2]=1;
    Point::order[4][3]=5;
    Point::order[4][4]=0;
    Point::order[4][5]=3;
    Point::order[4][6]=1;
    //
    Point::order[5][0]=1;
    Point::order[5][1]=2;
    Point::order[5][2]=0;
    Point::order[5][3]=5;
    Point::order[5][4]=1;
    Point::order[5][5]=3;
    Point::order[5][6]=0;
    //
    Point::order[6][0]=0;
    Point::order[6][1]=3;
    Point::order[6][2]=1;
    Point::order[6][3]=5;
    Point::order[6][4]=0;
    Point::order[6][5]=2;
    Point::order[6][6]=1;
    //
    Point::order[7][0]=1;
    Point::order[7][1]=3;
    Point::order[7][2]=0;
    Point::order[7][3]=5;
    Point::order[7][4]=1;
    Point::order[7][5]=2;
    Point::order[7][6]=0;
    //
    Point::phl.assign(5832,0);
    Point::dupes.assign(729,0);
    Point::corner_a.assign(27,0); //(0,0,1,0,0,1,2,2,3,0,0,1,0,0,1,2,2,3,4,4,5,4,4,5,6,6,7);
    Point::corner_b.assign(27,0); //(7,6,6,5,4,4,5,4,4,3,2,2,1,0,0,1,0,0,3,2,2,1,0,0,1,0,0);
    //
    Point::corner_a[0]=0;
    Point::corner_a[1]=0;
    Point::corner_a[2]=1;
    Point::corner_a[3]=0;
    Point::corner_a[4]=0;
    Point::corner_a[5]=1;
    Point::corner_a[6]=2;
    Point::corner_a[7]=2;
    Point::corner_a[8]=3;
    Point::corner_a[9]=0;
    Point::corner_a[10]=0;
    Point::corner_a[11]=1;
    Point::corner_a[12]=0;
    Point::corner_a[13]=0;
    Point::corner_a[14]=1;
    Point::corner_a[15]=2;
    Point::corner_a[16]=2;
    Point::corner_a[17]=3;
    Point::corner_a[18]=4;
    Point::corner_a[19]=4;
    Point::corner_a[20]=5;
    Point::corner_a[21]=4;
    Point::corner_a[22]=4;
    Point::corner_a[23]=5;
    Point::corner_a[24]=6;
    Point::corner_a[25]=6;
    Point::corner_a[26]=7;
    //
    Point::corner_b[0]=7;
    Point::corner_b[1]=6;
    Point::corner_b[2]=6;
    Point::corner_b[3]=5;
    Point::corner_b[4]=4;
    Point::corner_b[5]=4;
    Point::corner_b[6]=5;
    Point::corner_b[7]=4;
    Point::corner_b[8]=4;
    Point::corner_b[9]=3;
    Point::corner_b[10]=2;
    Point::corner_b[11]=2;
    Point::corner_b[12]=1;
    Point::corner_b[13]=0;
    Point::corner_b[14]=0;
    Point::corner_b[15]=1;
    Point::corner_b[16]=0;
    Point::corner_b[17]=0;
    Point::corner_b[18]=3;
    Point::corner_b[19]=2;
    Point::corner_b[20]=2;
    Point::corner_b[21]=1;
    Point::corner_b[22]=0;
    Point::corner_b[23]=0;
    Point::corner_b[24]=1;
    Point::corner_b[25]=0;
    Point::corner_b[26]=0;
      //
    Point::sequence.resize(8);
    Point::sequence[0]=0;      
    Point::sequence[1]=1;      
    Point::sequence[2]=3;      
    Point::sequence[3]=2;      
    Point::sequence[4]=6;      
    Point::sequence[5]=7;      
    Point::sequence[6]=5;      
    Point::sequence[7]=4;      
      //
    Point::updown.resize(8);
    Point::updown[0]=1;      
    Point::updown[1]=3;      
    Point::updown[2]=0;      
    Point::updown[3]=5;      
    Point::updown[4]=1;      
    Point::updown[5]=2;      
    Point::updown[6]=0;      
    Point::updown[7]=4;      
      //
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
	    //	    FileFractal << "candid " << pl << " " << pla << "   " << x << " " << y << " " << z << "  ";
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
			    //			    FileFractal << ph << " " << plb << " " << which+Point::dupes[pl*27+pla] << " " ;
			    Point::dupes[pl*27+pla]++;
			  }
		      }
		  }
	      }
	    //	    FileFractal << "\n";
	  }
      }
    int positions[27][27];
    bool decisions[27][27];
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
