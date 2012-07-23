#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void groups_level(Fractal& fractal,vector < vector<Group*> >& all_groups)
    {
      static ofstream GroupsPoints;
      if(!GroupsPoints.is_open())
	GroupsPoints.open("gp_lev.d");      
      vector <int> parts;
      parts.assign(100,0);
      vector <int> ngroups;
      vector <int> npoints;
      ngroups.assign(100,0);
      npoints.assign(100,0);
      double log1p5=log(1.5);
      int binmax=0;
      for(int ni=0;ni<fractal.get_number_particles();ni++)
	{
	  int lev=fractal.particle_list[ni]->get_highest_level();
	  parts[lev]++;
	}

      GroupsPoints << " " << endl;
      GroupsPoints << " steps " << fractal.get_steps() << endl;
      for(int ni=0;ni<=fractal.get_level_max();ni++)
	{
	  unsigned int sgroups=all_groups[ni].size();
	  unsigned int spoints=0;
	  for(unsigned int ng=0;ng<sgroups;ng++)
	    {
	      int ps=all_groups[ni][ng]->list_points.size();
	      spoints+=ps;
	      if(ni == 0) continue;
	      int bin=(double)(log((double)ps/26.999)/log1p5);
	      binmax=max(binmax,bin);
	      ngroups[bin]++;
	      npoints[bin]+=ps;
	    }
	  if(sgroups > 0) GroupsPoints << "\t" << ni << "\t" << sgroups << "\t" << spoints << "\t" << parts[ni] << endl;
	}
      GroupsPoints << endl;
      for(int ni=0;ni<=binmax;ni++)
	GroupsPoints << "\t" << ni << "\t" << (int)(27.0*pow(1.5,ni)+0.001) << "\t" << ngroups[ni] << 
	  "\t" << npoints[ni] << endl;
    }
}
