#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  typedef deque<double>::iterator _ITD__;
  template <class ForwardIterator>
  void start_writing(Fractal_Memory* PFM,int Numberparticles,
		     double G,vector <double>& xmin,vector <double>& xmax,
		     ForwardIterator posxb,ForwardIterator posyb,ForwardIterator poszb,
		     ForwardIterator velxb,ForwardIterator velyb,ForwardIterator velzb,ForwardIterator massesb)
  {
    static bool _DOIT=true;
    double t1=-PFM->p_mess->Clock();
    FILE* PFPos=PFM->p_file->PFPos;
    vector <double> pf(4);
    int RANK=PFM->p_mess->FractalRank;
    vector <int>BOX=PFM->Boxes[RANK];
    vector <double>RBOX=PFM->RealBoxes[RANK];
    const double conv_pot=G/(xmax[0]-xmin[0]);
    const double conv_force=conv_pot/(xmax[0]-xmin[0]);
    bool period=PFM->periodic;
    double timevar=PFM->time;
    if(period)
      timevar=PFM->arad;
    const double x0=-110.0;
    const double y0=120.0;
    const double z0=55.5;
    const double totalM=1.0e9;
    const double centerM=0.0;
    const double rmaX=30.0;
    
    const double slopE=-0.3;
    const double slopE2=slopE+2.0;
    const double slopE3=slopE+3.0;
    bool isoT=abs(slopE+2.0) < 0.01;
    const double consT=G*totalM/(pow(rmaX,slopE3)*slopE2);

    for(int ni=0;ni<Numberparticles;ni++)
      {
	PFM->p_fractal->particle_list[ni]->get_field_pf(pf);
	const double pot=pf[0]*conv_pot;
	const double fx=pf[1]*conv_force;
	const double fy=pf[2]*conv_force;
	const double fz=pf[3]*conv_force;
	int lev=PFM->p_fractal->particle_list[ni]->get_highest_level();
	fprintf(PFPos," Out %d %10.2E %6d L %d %13.6E %13.6E %13.6E",PFM->steps,timevar,ni,lev,*posxb,*posyb,*poszb); // 1-9
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E",*velxb,*velyb,*velzb,*massesb); // 10-13
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",pot,fx,fy,fz); // 14-17
 	const double dx=*posxb-x0;
 	const double dy=*posyb-y0;
 	const double dz=*poszb-z0;
 	const double dr=sqrt(dx*dx+dy*dy+dz*dz)+1.0e-50;
 	const double dr2=dr*dr;
 	const double dr3=dr2*dr;
 	const double fr=(dx*fx+dy*fy+dz*fz)/dr;
	const double ft=sqrt(pow(fr*dx/dr-fx,2)+pow(fr*dy/dr-fy,2)+pow(fr*dz/dr-fz,2));
	double mass=totalM+centerM;
	double potT=-G*centerM/dr;
	double fxT=-dx*G/dr3*centerM;
	double fyT=-dy*G/dr3*centerM;
	double fzT=-dz*G/dr3*centerM;
	double frTheory=-G*centerM/dr2;
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
	    potT+=-G*totalM/dr;
	  }
	fxT+=-dx*G/dr3*mass;
	fyT+=-dy*G/dr3*mass;
	fzT+=-dz*G/dr3*mass;
	if(_DOIT)
	  fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",potT,fxT,fyT,fzT); //18-21
	frTheory-=G*mass/dr2;
	const double difx=fx-fxT;
	const double dify=fy-fyT;
	const double difz=fz-fzT;
	const double ferror=sqrt(pow(difx,2)+pow(dify,2)+pow(difz,2));
	// const double err=ferror/abs(frTheory);
	const double err=abs((pot-potT)/potT);
	if(_DOIT)
	  fprintf(PFPos," EP %d %13.6E %13.6E %13.6E %13.6E %13.6E ",err > 0.1,dr,fr,frTheory,ft,ferror); // 22-28
	fprintf(PFPos," %13.6E %13.6E %13.6E ",dr,fr,ft); // 29-31
	if(_DOIT)
	  fprintf(PFPos," EF %d",abs((fr-frTheory)/frTheory) >0.1); // 32-33
	fprintf(PFPos,"\n");
	posxb++;
	posyb++;
	poszb++;
	velxb++;
	velyb++;
	velzb++;
      }
    // fflush(PFPos);
    t1+=PFM->p_mess->Clock();
    //    PFM->p_file->FileTime << " output time " << PFM->steps << " " << fixed << t1 << "\n";
    fprintf(PFM->p_file->PFTime," output time %5d %10.2E \n",PFM->steps,t1);
    massesb++;
    _DOIT=false;
  }
}
namespace FractalSpace
{
  template
  void start_writing(Fractal_Memory* PFM,int Numberparticles,
		     double G,vector <double>& xmin,vector <double>& xmax,
		     _ITD__ posxb,_ITD__ posyb,_ITD__ poszb,
		     _ITD__ velxb,_ITD__ velyb,_ITD__ velzb,_ITD__ massesb);
}
namespace FractalSpace
{
  void start_writing(Fractal_Memory* PFM,int Numberparticles,
		     double G,vector <double>& xmin,vector <double>& xmax,
		     vector<double>& posx,vector<double>& posy,vector<double>& posz,
		     vector<double>& velx,vector<double>& vely,vector<double>& velz,vector<double>& masses)
  {
    static bool _DOIT=true;
    double t1=-PFM->p_mess->Clock();
    FILE* PFPos=PFM->p_file->PFPos;
    vector <double> pf(4);
    int RANK=PFM->p_mess->FractalRank;
    vector <int>BOX=PFM->Boxes[RANK];
    vector <double>RBOX=PFM->RealBoxes[RANK];
    double conv_pot=G/(xmax[0]-xmin[0]);
    double conv_force=conv_pot/(xmax[0]-xmin[0]);
    bool period=PFM->periodic;
    double timevar=PFM->time;
    if(period)
      timevar=PFM->arad;
    double x0=-1.0;
    double y0=1.5;
    double z0=0.5;
    double totalM=1.0e9;
    double centerM=0.0;
    double rmaX=30.0;
    
    double slopE=-0.3;
    double slopE2=slopE+2.0;
    double slopE3=slopE+3.0;
    bool isoT=abs(slopE+2.0) < 0.01;
    double consT=G*totalM/(pow(rmaX,slopE3)*slopE2);

    for(int ni=0;ni<Numberparticles;ni++)
      {
	PFM->p_fractal->particle_list[ni]->get_field_pf(pf);
	double pot=pf[0]*conv_pot;
	double fx=pf[1]*conv_force;
	double fy=pf[2]*conv_force;
	double fz=pf[3]*conv_force;
	int lev=PFM->p_fractal->particle_list[ni]->get_highest_level();
	fprintf(PFPos," Out %d %10.2E %6d L %d %13.6E %13.6E %13.6E",PFM->steps,timevar,ni,lev,posx[ni],posy[ni],posz[ni]); // 1-9
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E",velx[ni],vely[ni],velz[ni],masses[ni]); // 10-13
	fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",pot,fx,fy,fz); // 14-17
 	double dx=posx[ni]-x0;
 	double dy=posy[ni]-y0;
 	double dz=posz[ni]-z0;
 	double dr=sqrt(dx*dx+dy*dy+dz*dz)+1.0e-50;
 	double dr2=dr*dr;
 	double dr3=dr2*dr;
 	double fr=(dx*fx+dy*fy+dz*fz)/dr;
	double ft=pow(fr*dx/dr-fx,2)+pow(fr*dy/dr-fy,2)+pow(fr*dz/dr-fz,2);
	ft=sqrt(ft);
	double mass=totalM+centerM;
	double potT=-G*centerM/dr;
	double fxT=-dx*G/dr3*centerM;
	double fyT=-dy*G/dr3*centerM;
	double fzT=-dz*G/dr3*centerM;
	double frTheory=-G*centerM/dr2;
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
	    potT+=-G*totalM/dr;
	  }
	fxT+=-dx*G/dr3*mass;
	fyT+=-dy*G/dr3*mass;
	fzT+=-dz*G/dr3*mass;
	if(_DOIT)
	  fprintf(PFPos," %13.6E %13.6E %13.6E %13.6E ",potT,fxT,fyT,fzT); //18-21
	frTheory-=G*mass/dr2;
	double difx=fx-fxT;
	double dify=fy-fyT;
	double difz=fz-fzT;
	double ferror=sqrt(pow(difx,2)+pow(dify,2)+pow(difz,2));
	// double err=ferror/abs(frTheory);
	double err=abs((pot-potT)/potT);
	if(_DOIT)
	  fprintf(PFPos," EP %d %13.6E %13.6E %13.6E %13.6E %13.6E ",err > 0.1,dr,fr,frTheory,ft,ferror); // 22-28
	fprintf(PFPos," %13.6E %13.6E %13.6E ",dr,fr,ft); // 29-31
	if(_DOIT)
	  fprintf(PFPos," EF %d",abs((fr-frTheory)/frTheory) >0.1); // 32-33
	fprintf(PFPos,"\n");
      }
    // fflush(PFPos);
    t1+=PFM->p_mess->Clock();
    //    PFM->p_file->FileTime << " output time " << PFM->steps << " " << fixed << t1 << "\n";
    fprintf(PFM->p_file->PFTime," output time %5d %10.2E \n",PFM->steps,t1);
    _DOIT=false;
  }
}
