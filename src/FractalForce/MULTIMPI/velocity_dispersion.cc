#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void velocity_dispersion(
			   vector <double>& radius,
			   vector <double>& vel_disp,
			   const int& datapoints,
			   const int& profile,
			   const double& xmin,
			   const double& xmax,
			   const double& slope,
			   const double& scale_factor,
			   const double& vel_factor)
  {
    radius.resize(datapoints);
    vel_disp.resize(datapoints);
    vector <double>rhosig2(datapoints,0.0);
    const int counts=1000;
    const double xlast=1.0e6*xmax;
  }
  double integrate_lin(const double& xa,
		       const double& xb,
		       const int& counts)
  {
    for(int ni=1;ni<counts;ni++)
      {
	
      }
  }
  double integrate_exp(const double& xa,
		       const double& xb,
		       const int& counts)
  {
    
  }
}
