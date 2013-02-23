#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void start_writing(Fractal_Memory* PFM,int Numberparticles,
		     vector<double>& posx,vector<double>& posy,vector<double>& posz,
		     vector<double>& velx,vector<double>& vely,vector<double>& velz,vector<double>& masses)
  {
    ofstream& FilePos=PFM->p_file->FilePos;
    int NP=PFM->number_particles;
    vector <double> pf(4);
    for(int ni=0;ni<NP;ni++)
      {
	PFM->p_fractal->particle_list[ni]->get_field_pf(pf);
	FilePos << "out " << PFM->steps << "\t" << ni << "\t" << fixed << posx[ni] << "\t" << posy[ni] << "\t" << posz[ni];
	FilePos << scientific << "\t" << velx[ni] << "\t" << vely[ni] << "\t" << velz[ni];
	FilePos << "\t" << masses[ni];
	FilePos << scientific << "\t" << pf[0] << "\t" << pf[1] << "\t" << pf[2] << "\t" << pf[3] ;
	FilePos  << endl;
      }
  }
}
