#include "iostream.h"
#include "math.h"
//
float Age_of_the_universe (float omega_0, float omega_lambda, float redshift);
float Omega (float omega_0, float omega_lambda, float redshift);
float Lambda (float omega_0, float omega_lambda, float redshift);
float Hubble (float omega_0, float omega_lambda, float redshift);
float Growth(float omega_0, float omega_lambda, float redshift);
float Length(float omega_0, float omega_lambda, float redshift);
//
  int main()
{
  float omega_0,omega_lambda,redshift;
  cout << "Omega_0, Omega_lambda, redshift " ;
  cin >> omega_0;
  cin  >> omega_lambda;
  cin >> redshift;
  cout << "\n";
  cout << "Omega_0 \t"  << omega_0 << "\n";
  cout << "Omega_lambda \t" << omega_lambda << "\n";
  cout << "Redshift \t" << redshift << "\n";
  cout << "Age(z)  \t" <<  Age_of_the_universe(omega_0, omega_lambda, redshift) << "\n";
  cout << "Hubble(z) \t" <<  Hubble(omega_0, omega_lambda,redshift) << "\n";
  cout << "Omega(z) \t" <<  Omega(omega_0, omega_lambda,redshift) << "\n";
  cout << "Lambda(z) \t" <<  Lambda(omega_0, omega_lambda,redshift) << "\n";
  cout << "Length(z) \t" <<  Length(omega_0, omega_lambda,redshift) << "\n";
  cout << "Growth(z)  \t" <<  Growth(omega_0, omega_lambda,redshift) << "\n";
  //
  float  growth_1=Growth(omega_0, omega_lambda,redshift+0.01);
  float  growth_2=Growth(omega_0, omega_lambda,redshift-0.01);
  cout << "Growth Rate(z) \t " << log(growth_2/growth_1)/log((1.0+redshift+0.01)/(1.0+redshift-0.01)) << "\n";
  //
  return 0;
}
//
//
float Hubble (float omega_0, float omega_lambda, float redshift)
{
//
  return sqrt(omega_0*pow(1.0+redshift,3)+
	      (1.0-omega_0-omega_lambda)*pow(1.0+redshift,2)+omega_lambda);
}
//
float Omega (float omega_0, float omega_lambda, float redshift)
//
{
  return omega_0*pow(1.0+redshift,3)/
    pow(Hubble(omega_0,omega_lambda,redshift),2);  
}
//
float Lambda(float omega_0, float omega_lambda, float redshift)
//
{
  return omega_lambda/pow(Hubble(omega_0,omega_lambda,redshift),2);  
}
//
float Age_of_the_universe(float omega_0, float omega_lambda, float redshift)
{
  float age=0.0;
  float dx;
  float x;
  int n;
  //
  dx=1.0e-4/(1.0+redshift);
  //
  for (n=0;n < 10000;n++)
    {
      x=dx*n+dx*0.5;
      age+=dx*x/sqrt(omega_0*x+(1.0-omega_0-omega_lambda)*x*x+
		     omega_lambda*pow(x,4));
    }
  //
  return age;
}
//
float Growth(float omega_0, float omega_lambda, float redshift)
//
{
  float omega= Omega (omega_0, omega_lambda, redshift);
  //
  float lambda= Lambda (omega_0, omega_lambda, redshift);
  //
  float g=2.5*omega/(pow(omega,4.0/7.0)-lambda+
		     (1.0+omega*0.5)*(1.0+lambda/70.0));
  //
  float g_0=2.5*omega_0/(pow(omega_0,4.0/7.0)-omega_lambda+
		     (1.0+omega_0*0.5)*(1.0+omega_lambda/70.0));
  return g/(1.0+redshift)/g_0;
}
//
float Length(float omega_0, float omega_lambda, float redshift)
{
  float dist=0.0;
  float dx;
  float x;
  float start;
  int n;
  //
  start=1.0/(1.0+redshift);
  dx=1.0e-4*(1.0-start);
  //
  for (n=0;n < 10000;n++)
    {
      x=dx*n+dx*0.5+start;
      dist+=dx/(x*x*Hubble (omega_0, omega_lambda, 1.0/x-1.0));
    }
  //
  if(omega_0 < 1.0)
    dist=sinh(sqrt(1.0-omega_0)*dist)/sinh(sqrt(1.0-omega_0));
  else if(omega_0 > 1.0)
    dist=sin(sqrt(omega_0-1.0)*dist)/sin(sqrt(omega_0-1.0));
  //
  return dist;
}
