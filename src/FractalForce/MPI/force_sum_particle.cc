#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_sum_particle(Group& group,Fractal& fractal,const bool& doit)
  { 
    //
    vector <double> pf(4);
    vector <double> sumf(3);
    sumf.assign(3,0.0);
    vector <double> sumf0(3);
    sumf0.assign(3,0.0);
    double sum0=1.0e-30;
    if(doit) 
      {
	group.get_forcem(sumf0,sum0);
	//	cout << " SUM0 " << &group << " " << sumf0[0] << " " << sumf0[1] << " " << sumf0[2] << " " << sum0 << endl;
	sum0=1.0e-30;
      }
    //
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(point.list_particles.empty()) continue;
	//
	for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
	  {
	    Particle& particle=**particle_itr;
	    if(particle.get_p_highest_level_group() == 0) continue;
	    particle.get_field_pf(pf);
	    double mass=particle.get_mass();
	    sumf[0]+=pf[1]*mass;
	    sumf[1]+=pf[2]*mass;
	    sumf[2]+=pf[3]*mass;
	    sum0+=mass;
	  }
      }
    if(doit)
      {
	pf[0]=0.0;
	pf[1]=(sumf0[0]-sumf[0])/sum0;
	pf[2]=(sumf0[1]-sumf[1])/sum0;
	pf[3]=(sumf0[2]-sumf[2])/sum0;
	//	cout << " SUMM " << &group << " " << sumf[0] << " " << sumf[1] << " " << sumf[2] << " " << sum0 << endl;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    Point& point=**point_itr;
	    if(point.list_particles.empty()) continue;
	    //
	    for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
	      {
		Particle& particle=**particle_itr;
		if(particle.get_p_highest_level_group() == 0) continue;
		particle.add_field_pf(pf);
	      }
	  }
      }	
    else
      {
	group.set_forcem(sumf,sum0);
	//	cout << " SUMS " << &group << " " << sumf[0] << " " << sumf[1] << " " << sumf[2] << " " << sum0 << endl;
      }
  }
}
