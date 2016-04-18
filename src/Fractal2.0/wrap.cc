#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void wrap(const int& p)
  {
    vector <double> pos(3);
    wrap(particle_list[p]);
  }
  void wrap()
  {
    for(unsigned int p=0;p<particle_list.size();p++)
      wrap(p);
  }
  void wrap(Particle* p)
  {
    vector <double> pos(3);
    p->get_pos(pos);
    bool doit=false;
    for(int j=0;j<3;j++)
      {
	if(pos[j] < 0.0)
	  {
	    pos[j]++;
	    doit=true;
	  }
	else if(pos[j] >= 1.0)
	  {
	    pos[j]--;
	    doit=true;
	  }
      }
    if(doit)
      p->set_pos(pos);
  }
}
