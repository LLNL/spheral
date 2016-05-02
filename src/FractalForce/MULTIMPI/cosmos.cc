
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  double Hubble (const double& omega_0, const double& omega_lambda, const double& redshift)
  {
    //
    return sqrt(omega_0*pow(1.0+redshift,3)+
		(1.0-omega_0-omega_lambda)*pow(1.0+redshift,2)+omega_lambda);
  }
  //
  double Omega (const double& omega_0, const double& omega_lambda, const double& redshift)
  //
  {
    return omega_0*pow(1.0+redshift,3)/
      pow(Hubble(omega_0,omega_lambda,redshift),2);  
  }
  //
  double Lambda(const double& omega_0, const double& omega_lambda, const double& redshift)
  //
  {
    return omega_lambda/pow(Hubble(omega_0,omega_lambda,redshift),2);  
  }
  //
  double Age_of_the_universe(const double& omega_0, const double& omega_lambda, const double& redshift)
  {
    double age=0.0;
    double dx=1.0e-4/(1.0+redshift);
    //
    for (int n=0;n < 10000;n++)
      {
	double x=dx*n+dx*0.5;
	age+=dx*x/sqrt(omega_0*x+(1.0-omega_0-omega_lambda)*x*x+
		       omega_lambda*pow(x,4));
      }
    //
    return age;
  }
  //
  double Growth(const double& omega_0, const double& omega_lambda, const double& redshift)
  //
  {
    double omega= Omega (omega_0, omega_lambda, redshift);
    //
    double lambda= Lambda (omega_0, omega_lambda, redshift);
    //
    double g=2.5*omega/(pow(omega,4.0/7.0)-lambda+
			(1.0+omega*0.5)*(1.0+lambda/70.0));
    //
    double g_0=2.5*omega_0/(pow(omega_0,4.0/7.0)-omega_lambda+
			    (1.0+omega_0*0.5)*(1.0+omega_lambda/70.0));
    return g/(1.0+redshift)/g_0;
  }
  double dGrowthdT(const double& omega_0, const double& omega_lambda, const double& redshift)
  //
  {
    double dGdA=dGrowthdA(omega_0, omega_lambda,redshift);
    double dAdT=Hubble(omega_0,omega_lambda,redshift)/(1.0+redshift);
    return dGdA*dAdT;
  }
  double dGrowthdA(const double& omega_0, const double& omega_lambda, const double& redshift)
  //
  {
    double arad=1.0/(redshift+1.0);
    double aradup=arad+0.01;
    double araddown=arad-0.01;
    double dup=Growth(omega_0,omega_lambda,1.0/aradup-1.0);
    double ddown=Growth(omega_0,omega_lambda,1.0/araddown-1.0);
    return (dup-ddown)/0.02;
  }
  //
  double Length(const double& omega_0, const double& omega_lambda, const double& redshift)
  {
    double dist=0.0;
    double start=1.0/(1.0+redshift);
    double dx=1.0e-4*(1.0-start);
    //
    for (int n=0;n < 10000;n++)
      {
	double x=dx*n+dx*0.5+start;
	dist+=dx/(x*x*Hubble (omega_0, omega_lambda, 1.0/x-1.0));
      }
    //
    if(omega_0 < 0.9999)
      dist=sinh(sqrt(1.0-omega_0)*dist)/sinh(sqrt(1.0-omega_0));
    else if(omega_0 > 1.0001)
      dist=sin(sqrt(omega_0-1.0)*dist)/sin(sqrt(omega_0-1.0));
    //
    return dist;
  }
}
