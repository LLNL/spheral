#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void tree_dump(Fractal_Memory& FM)
  {
    //    ofstream& FT=FM.p_file->FileTree;
    //    ofstream& FS=FM.p_file->FileSurface;
    FILE* PFT=FM.p_file->PFTree;
    FILE* PFS=FM.p_file->PFSurface;
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
		    //		    FT << " S " << FM.steps << " I " << inside << " L " << level << " G " << gnumber << " " << pnumber << " " << dens;
		    //		    FT << " " << pos[0]  << " " << pos[1]  << " " << pos[2] << "\n";
		    fprintf(PFT,"S%d I%d L%2d G%5d %6d %13.6E",FM.steps,inside,level,gnumber,pnumber,dens);
		    fprintf(PFT,"%8d %8d %8d \n",pos[0],pos[1],pos[2]);
		  }
		else if(level > 0)
		  {
		    for(int ni=1;ni<6;ni+=2)
		      {
			Point* p=point.get_point_ud(ni);
			if(p && !p->get_inside())
			  {
			    p->get_pos_point(posup);
			    fprintf(PFS,"S%d L%2d G%5d %6d",FM.steps,level,gnumber,pnumber);
			    fprintf(PFS,"%8d %8d %8d",pos[0],pos[1],pos[2]);
			    fprintf(PFS,"%8d %8d %8d \n",posup[0],posup[1],posup[2]);
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
