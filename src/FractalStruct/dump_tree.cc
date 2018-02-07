#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void dump_tree(Fractal_Memory& fractal_memory,Fractal& fractal)
  {
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    FileFractal << "enter dump_tree " << "\n";
    vector <int>Box(6);
    fractal.getBox(Box);
    int zoom=Misc::pow(2,fractal.get_level_max());
    for(int i=0;i<6;i++)
      Box[i]=Box[i]*zoom;
    for(int lev=0;lev <= fractal.get_level_max();lev++)
      {
	FileFractal << "dump level " << lev << "\n";
	// for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[lev].begin();
	//     group_itr!=fractal_memory.all_groups[lev].end();group_itr++)
	for(auto pg : fractal_memory.all_groups[lev])
	  {
	    Group& group=*pg;
	    // Group& group=**group_itr;
	    FileFractal << "dump group " << &group << "\n";
	    // for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
	    for(auto pp : group.list_points)
	      {
		Point& point=*pp;
		// Point& point=**point_itr;
		vector <int> ppos(3);
		point.get_pos_point(ppos);
		bool inside=      
		  ppos[0] >= Box[0] && ppos[0] <= Box[1] &&
		  ppos[1] >= Box[2] && ppos[1] <= Box[3] &&
		  ppos[2] >= Box[4] && ppos[2] <= Box[5];
		bool pop=point.list_particles.size() > 0;
		if(inside == pop) continue;
		point.dump();
		// for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
		for(auto part : point.list_particles)
		  // {
		  //   Particle& particle=**particle_itr;
		    // particle.dump(FileFractal);
		  part->dump(FileFractal);
		  // }
	      }
	  }
      }
  }
}
