#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void balance_by_points(Fractal_Memory& mem,Fractal& frac)
  {
    int FractalNodes0=mem.p_mess->FractalNodes0;
    int FractalNodes1=mem.p_mess->FractalNodes1;
    int FractalNodes2=mem.p_mess->FractalNodes2;
    int Root=0;
    int length=32;
    int length3=length*length*length;
    double alength=length;
    int* numbers=new int[length3];
    for(int ni=0;ni<length3;ni++)
      numbers[ni]=0;
    vector <double> pos(3);
    for(int n=0;n < frac.get_number_particles();++n)
      {
	Particle* p=frac.particle_list[n];
	p->get_pos(pos);
	double ax=pos[0]*alength;
	if(ax < 0.0 || ax >= alength) continue;
	double ay=pos[1]*alength;
	if(ay < 0.0 || ay >= alength) continue;
	double az=pos[2]*alength;
	if(az < 0.0 || az >= alength) continue;
	int nx=ax;
	int ny=ay;
	int naz=az;
	int n=nx+(ny+nz*length)*length;
	numbers[n]++;
      }
    mem.p_mess->Find_Sum_INT_to_Root(numbers,length3,Root);
    mem.p_mess->Send_INT_from_Root(numbers,length3,Root);
    vector <int>points(length3)
    for(int ni=0;ni<length3;ni++)
      {
	numbers[ni]=min(numbers[ni],1);
	points[ni]=numbers[ni];
      }
    vector <int>pointz(length+1);
    pointz[0]=0;
    for(int nz=0;nz<length;nz++)
      {
	pointz[nz+1]=pointz[nz];
	for(int ny=0;ny<length;ny+)
	  {
	    for(int nx=0;nx<length;nx+)
	      {
		int n=nx+(ny+nz*length)*length;
		pointz[nz+1]+=points[ni];
	      }
	  }
      }
    int pztotal=pointz[length];
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	int target=(FRZ*pztotal)/FractalNodes2;
	for(int nz=0;nz<length;nz++)
	  {
	    if(pointz[nz] >= target)
	      {
		success=true;
		boxz[FRZ]=
		
	      }
	  }
	}
    delete [] numbers;
  }
}
