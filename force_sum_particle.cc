#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_sum_particle(Group& group,const bool& doit)
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
	//	cerr << " SUM0 " << &group << " " << sumf0[0] << " " << sumf0[1] << " " << sumf0[2] << " " << sum0 << "\n";
	sum0=1.0e-30;
      }
    //
    for(auto &p : group.list_points)
      {
	if(p->list_particles.empty())
	  continue;
	//
	for(auto &part : p->list_particles)
	  {
	    if(part->get_p_highest_level_group() == 0)
	      continue;
	    part->get_field_pf(pf);
	    double mass=part->get_mass();
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
	//	cerr << " SUMM " << &group << " " << sumf[0] << " " << sumf[1] << " " << sumf[2] << " " << sum0 << "\n";
	for(auto pp : group.list_points)
	  for(auto part : pp->list_particles)
	    if(part->get_p_highest_level_group())
	      part->add_field_pf(pf);
      }	
    else
      {
	group.set_forcem(sumf,sum0);
	//	cerr << " SUMS " << &group << " " << sumf[0] << " " << sumf[1] << " " << sumf[2] << " " << sum0 << "\n";
      }
  }
}
