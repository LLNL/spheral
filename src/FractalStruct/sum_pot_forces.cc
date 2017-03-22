#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void sum_pot_forces(Fractal& frac)
  {
    ofstream& FileFractal=frac.p_file->DUMPS;
    //    ofstream& FileFractal=frac.p_file->FileFractal;
    double n=frac.get_number_particles();
    vector <double> field(4);
    vector <double> sum_8(8,0.0);
    for(int p=0;p < n;++p)
      {
	if(!frac.particle_list[p]->get_real_particle())
	  continue;
	frac.particle_list[p]->get_field_pf(field);
	for(int ni=0;ni<4;ni++)
	  {
	    sum_8[ni*2]+=field[ni];
	    sum_8[ni*2+1]+=field[ni]*field[ni];
	  }
      }    
    for(int ni=0;ni < 4;ni++)
      {
	sum_8[ni*2]/=n;
	sum_8[ni*2+1]/=n;
	sum_8[ni*2+1]=sqrt(sum_8[ni*2+1]-sum_8[ni*2]*sum_8[ni*2]);
      }
    FileFractal << "sum_start " << sum_8[0] << " " << sum_8[1] << " " << sum_8[2] << " " << sum_8[3] << " " << sum_8[4] << " " << sum_8[5] << " " << sum_8[6] << " " << sum_8[7] << "\n";
  }
}
