#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void test_points(Fractal_Memory& mem,vector<vector<Point*>>& SPoints,int level)
  {
    static int levelmax=-1;
    if(level < levelmax)
      {
	levelmax=100;
	return;
      }
    levelmax=level;
    vector<int>BOX=mem.BoxesLev[mem.p_mess->FractalRank][level];
    ofstream& FHT=mem.p_file->DUMPS;
    double zoom=(double)mem.p_fractal->get_grid_length()*Misc::pow(2,mem.p_fractal->get_level_max());
    vector <double> xmin(3,-50.0);
    vector <double> xmax(3,50.0);
    double scalex=(xmax[0]-xmin[0])/zoom;
    double G=3.141592;
    double x0=-1.0;
    double y0=1.5;
    double z0=0.5;
    double totalM=1.0e9;
    double rmaX=30.0;
    double slopE=-1.8;
    double slopE2=slopE+2.0;
    double slopE3=slopE+3.0;
    // bool isoT=abs(slopE+2.0) < 0.01;
    double consT=G*totalM/(pow(rmaX,slopE3)*slopE2);
    double scalinga=G/(xmax[0]-xmin[0]);
    // for(auto pgroup : mem.all_groups[level])
    //   {
    // 	for(auto p : pgroup->list_points)
    // 	  {
    // 	    if(p->get_passive_point())
    // 	      continue;
    for(auto SP : SPoints)
      {
	for(auto p : SP)
	  {
	    if(!p)
	      continue;
	    array<int,3>pos=p->get_pos_point_a();
	    double pot=scalinga*p->get_potential_point();
	    double x=(double)pos[0]*scalex+xmin[0];
	    double y=(double)pos[1]*scalex+xmin[1];
	    double z=(double)pos[2]*scalex+xmin[2];
	    double dx=x-x0;
	    double dy=y-y0;
	    double dz=z-z0;
	    double dr=sqrt(dx*dx+dy*dy+dz*dz);
	    double potT=consT*(pow(dr,slopE2)-pow(rmaX,slopE2));
	    potT-=G*totalM/rmaX;
	    double eror=abs(potT-pot)/abs(potT);
	    bool onedge=on_edge(pos,BOX);
	    bool inside=p->get_inside();
	    bool trouble=p->get_trouble();
	    FHT << " POINTS " << level << " " << dr << " " << pot << " " << potT << " " << abs(pot) << " " << abs(potT) << " " << eror << " " << onedge << " " << inside << " " << trouble;
	    if(eror > 0.05)
	      {
		FHT << " EROR " << pos[0] << " " << pos[1] << " " << pos[2] << " ";
		FHT << BOX[0] << " " << BOX[1] << " " << BOX[2] << " " << BOX[3] << " " << BOX[4] << " " << BOX[5];
	      }
	
	      FHT << "\n";
	  }
      }
  }
}
