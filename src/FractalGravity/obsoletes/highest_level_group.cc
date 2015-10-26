#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void highest_level_group(list <Chain*>& list_chains,Fractal& fractal,Misc& misc)
  {
    cout << " enter highest " << endl;
    for(list <Chain*>::const_iterator chain_itr=list_chains.begin();chain_itr != list_chains.end();++chain_itr)
      {
	Chain& chain=**chain_itr;
	//	cout << "chain " << &chain << endl;
	for(list <Group*>::const_iterator group_itr=chain.list_groups.begin();group_itr!=chain.list_groups.end();++group_itr)
	  {
	    Group& group=**group_itr;
	    //	    cout << "group h " << &group << " " << misc.p_group_0 << " " << group.get_level() << endl;
	    if(&group != misc.p_group_0)
	      {
		int lev=group.get_level();
		assert(lev > 0);
		for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr!=group.list_points.end();++point_itr)
		  {
		    Point& point=**point_itr;
		    //		    cout << "point " << &point << endl;
		    for(list <Particle*>::const_iterator p_itr=point.list_particles.begin();p_itr!=point.list_particles.end();++p_itr)
		      {
			Particle& particle=**p_itr;
			//			cout << "particle " << particle << endl;
			int highest_level=particle.get_highest_level();
			assert(highest_level != lev);
			if(lev > highest_level)
			  {
			    particle.set_p_highest_level_group(&group);
			    particle.set_highest_level(lev);
			  }
		      }
		  }
	      }
	  }
      }
    cout << " exit highest " << endl;
  }
}
