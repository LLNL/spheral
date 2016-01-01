#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void start_writing(Fractal_Memory* PFM,int Numberparticles,
		     double G,vector <double>& xmin,vector <double>& xmax,
		     vector<double>& posx,vector<double>& posy,vector<double>& posz,
		     vector<double>& velx,vector<double>& vely,vector<double>& velz,vector<double>& masses)
  {
    double t1=-PFM->p_mess->Clock();
    FILE* PFPos=PFM->p_file->PFPos;
    vector <double> pf(4);
    double conv_pot=G/(xmax[0]-xmin[0]);
    double conv_force=conv_pot/(xmax[0]-xmin[0]);
    bool period=PFM->periodic;
    double timevar=PFM->time;
    if(period)
      timevar=PFM->arad;

    //    double x0=1.0;
    //    double y0=-2.0;
    //    double z0=3.0;
    //    double totalM=1.0e9;
    //    double centerM=0.0;
    //    double rmaX=30.0;

    //    double slopE=-1.5;
    //    double slopE2=slopE+2.0;
    //    double slopE3=slopE+3.0;
    //    bool isoT=abs(slopE+2.0) < 0.01;
    //    double consT=G*totalM/(pow(rmaX,slopE3)*slopE2);

    for(int ni=0;ni<Numberparticles;ni++)
      {
	PFM->p_fractal->particle_list[ni]->get_field_pf(pf);
	int lev=PFM->p_fractal->particle_list[ni]->get_highest_level();
	fprintf(PFPos," Out%d %10.2E %6d L%d %13.6E %13.6E %13.6E",PFM->steps,timevar,ni,lev,posx[ni],posy[ni],posz[ni]);
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E",velx[ni],vely[ni],velz[ni],masses[ni]);
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",pf[0]*conv_pot,pf[1]*conv_force,pf[2]*conv_force,pf[3]*conv_force);
	/*
	double dx=posx[ni]-x0;
	double dy=posy[ni]-y0;
	double dz=posz[ni]-z0;
	double dr=sqrt(dx*dx+dy*dy+dz*dz)+1.0e-50;
	double dr2=dr*dr;
	double dr3=dr2*dr;
	double fr=-conv_force*(dx*pf[1]+dy*pf[2]+dz*pf[3])/dr;
	double mass=totalM+centerM;
	double potT=-G*centerM/dr;
	double fxT=-dx*G/dr3*centerM;
	double fyT=-dy*G/dr3*centerM;
	double fzT=-dz*G/dr3*centerM;
	double frTheory=G*centerM/dr2;
	if(dr < rmaX)
	  {
	    mass=pow(dr/rmaX,slopE3)*totalM;
	    if(isoT)
	      potT+=totalM/rmaX*log(dr/rmaX);
	    else
	      potT+=consT*(pow(dr,slopE2)-pow(rmaX,slopE2));
	    potT-=G*totalM/rmaX;
	  }
	else
	  {
	    potT-=-G*totalM/dr;
	  }
	fxT+=-dx*G/dr3*mass;
	fyT+=-dy*G/dr3*mass;
	fzT+=-dz*G/dr3*mass;
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",potT,fxT,fyT,fzT);
	frTheory+=G*mass/dr2;
	double ferror=sqrt(pow(fxT-pf[1]*conv_force,2)+pow(fyT-pf[2]*conv_force,2)+pow(fzT-pf[3]*conv_force,2));
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",dr,fr,frTheory,ferror);
	*/
	fprintf(PFPos,"\n");
      }
    //    fflush(PFPos);
    t1+=PFM->p_mess->Clock();
    //    PFM->p_file->FileTime << " output time " << PFM->steps << " " << fixed << t1 << "\n";
    fprintf(PFM->p_file->PFTime," output time %5d %10.2E \n",PFM->steps,t1);
  }
}
