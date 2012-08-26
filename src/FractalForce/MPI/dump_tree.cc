#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void dump_tree(Fractal_Memory& fractal_memory,Fractal& fractal)
  {
    cout << "enter dump_tree " << endl;
    vector <int>Box(6);
    fractal.getBox(Box);
    int zoom=Misc::pow(2,fractal.get_level_max());
    for(int i=0;i<6;i++)
      Box[i]=Box[i]*zoom;
    for(int lev=0;lev <= fractal.get_level_max();lev++)
      {
	cout << "dump level " << lev << endl;
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[lev].begin();
	    group_itr!=fractal_memory.all_groups[lev].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    cout << "dump group " << &group << endl;
	    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
	      {
		Point& point=**point_itr;
		vector <int> ppos(3);
		point.get_pos_point(ppos);
		bool inside=      
		  ppos[0] >= Box[0] && ppos[0] <= Box[1] &&
		  ppos[1] >= Box[2] && ppos[1] <= Box[3] &&
		  ppos[2] >= Box[4] && ppos[2] <= Box[5];
		bool pop=point.list_particles.size() > 0;
		if(inside == pop) continue;
		point.dump();
		for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
		  {
		    Particle& particle=**particle_itr;
		    particle.dump();
		  }
	      }
	  }
      }
  }
}
