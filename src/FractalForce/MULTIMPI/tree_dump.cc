#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void tree_dump(Fractal_Memory& FM)
  {
    ofstream& FT=FM.p_file->FileTree;
    ofstream& FS=FM.p_file->FileSurface;
    vector <int>pos(3);
    vector <int>posup(3);
    for(int level=0;level <= FM.level_max;level++)
      {
	int gnumber=0;
	for(vector <Group*>::const_iterator group_itr=FM.all_groups[level].begin();
	    group_itr!=FM.all_groups[level].end();group_itr++)
	  {
	    int pnumber=0;
	    Group& group=**group_itr;
	    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point& point=**point_itr;
		point.get_pos_point(pos);
		double dens=point.get_density_point();
		bool inside=point.get_inside();
		if(inside)
		  {
		    FT << " S " << FM.steps << " I " << inside << " L " << level << " G " << gnumber << " " << dens;
		    FT << " " << pos[0]  << " " << pos[1]  << " " << pos[2] << endl;
		  }
		else if(level > 0)
		  {
		    for(int ni=1;ni<6;ni+=2)
		      {
			Point* p=point.get_point_ud(ni);
			if(p && !p->get_inside())
			  {
			    p->get_pos_point(posup);
			    FS << " S " << FM.steps << " L " << level << " G " << gnumber;
			    FS << " " << pos[0]  << " " << pos[1]  << " " << pos[2];
			    FS << " " << posup[0] << " " << posup[1] << " " << posup[2] << endl;
			  }
		      }
		  }
		pnumber++;
	      }
	    gnumber++;
	  }
      }
  }
}
