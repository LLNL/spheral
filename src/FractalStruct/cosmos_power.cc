#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::power_spec="cdm";
  //
  double cosmos_power(const double& q,Fractal_Memory& fractal_memory)
  {
    //    ofstream& FileFractal=mem.p_fractal->p_file->FileFractal;
    double amplitude=pow(q,fractal_memory.power_slope)*exp(-pow(q/fractal_memory.cut_off,2));
    //    FileFractal << "cosmo " << q << " " << fractal_memory.cut_off << " " << fractal_memory.power_slope << " " << amplitude << "\n";
    if(fractal_memory.spectrum_number == 0)
      {
	return amplitude;
      }
    else if(fractal_memory.spectrum_number ==1)
      {
	double t=log(1.0+2.34*q)*pow(1.0+13.0*q+pow(10.5*q,2)+pow(10.4*q,3)+pow(6.51*q,4),-0.25)/(2.34*q);
	return amplitude*t*t;
      }
    else
      {
	assert(0);
      }
    return 0.0;
  }
}
