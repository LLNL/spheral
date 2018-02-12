#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void tree_dump(Fractal_Memory& FM)
  {
    //    FILE* PFT=FM.p_file->PFTree;
    FILE* PFS=FM.p_file->PFSurface;
    vector <int>pos(3);
    vector <int>posup(3);
    vector <Point*>ud;
    for(int level=1;level <= FM.level_max;level++)
      {
	int gnumber=0;
	// for(vector <Group*>::const_iterator group_itr=FM.all_groups[level].begin();
	//     group_itr!=FM.all_groups[level].end();group_itr++)
	for(auto pg : FM.all_groups[level])
	  {
	    int pnumber=0;
	    // Group& group=**group_itr;
	    // for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	    for(auto pp : pg->list_points)
	      {
		// Point& point=**point_itr;
		if(!pp->get_inside() && !pp->get_passive_point())
		  {
		    pp->get_pos_point(pos);
		    pp->get_point_ud(ud);
		    for(int ni=1;ni<6;ni+=2)
		      {
			if(ud[ni]==0 || ud[ni]->get_inside() || ud[ni]->get_passive_point())
			  continue;
			ud[ni]->get_pos_point(posup);
			fprintf(PFS,"SStep%5d L%3d G%7d",FM.steps,level,gnumber);
			fprintf(PFS,"%8d %8d %8d %8d %8d %8d \n",pos[0],pos[1],pos[2],posup[0]-pos[0],posup[1]-pos[1],posup[2]-pos[2]);
		      }
		  }
		pnumber++;
	      }
	    gnumber++;
	  }
      }
  }
}
