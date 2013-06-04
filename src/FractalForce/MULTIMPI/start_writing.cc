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
    ofstream& FilePos=PFM->p_file->FilePos;
    int NP=PFM->number_particles;
    vector <double> pf(4);
    double conv_pot=G/(xmax[0]-xmin[0]);
    double conv_force=conv_pot/(xmax[0]-xmin[0]);
    bool period=PFM->periodic;
    double timevar=PFM->time;
    if(period)
      timevar=PFM->arad;
    for(int ni=0;ni<NP;ni++)
      {
	PFM->p_fractal->particle_list[ni]->get_field_pf(pf);
	int lev=PFM->p_fractal->particle_list[ni]->get_highest_level();
	FilePos << "out " << PFM->steps << "\t" << timevar << "\t" << ni << "\t" << "L" << lev << "\t" << scientific << posx[ni] << "\t" << posy[ni] << "\t" << posz[ni];
	FilePos << scientific << "\t" << velx[ni] << "\t" << vely[ni] << "\t" << velz[ni];
	FilePos << "\t" << masses[ni];
	FilePos << scientific << "\t" << pf[0]*conv_pot << "\t" << pf[1]*conv_force << "\t" << pf[2]*conv_force << "\t" << pf[3]*conv_force ;
	FilePos  << endl;
      }
  }
}
