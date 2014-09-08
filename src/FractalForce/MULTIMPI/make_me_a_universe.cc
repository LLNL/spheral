#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_me_a_universe(Fractal_Memory* PFM,vector <double>& xmin,vector <double>& xmax,
			  vector <double>& posx,vector <double>& posy,vector <double>& posz,
			  vector <double>& velx,vector <double>& vely,vector <double>& velz,
			  vector <double>& masses)
  {
    Fractal* PF=PFM->p_fractal;
    ofstream& FF=PFM->p_file->FF;
    int FractalNodes=PFM->p_mess->FractalNodes;
    int FractalRank=PFM->p_mess->FractalRank;
    std::srand(9973+256*FractalRank);
    xmin.assign(3,0.0);
    xmax.assign(3,1.0);
    int total_number_particles=mem.p_mess->total_number_particles;
    double pi=4.0*atan(1.0);
    double total_mass=0.375/pi*PFM->omega_start;
    double G=1.0;

    double m=0.375*PFM->omega_start/pi/static_cast<double>(PFM->p_mess->number_particles_total);
    FF << "m= " << m << "\n";
    bool zel_predict=PFM->crash_levels > 0 && PFM->max_particles > PFM->number_particles; 
    int splits_tmp=0;
    double cut_off_tmp=0.0;
    if(zel_predict)
      {
	splits_tmp=PFM->splits;
	PFM->splits=0;
	cut_off_tmp=PFM->cut_off;
	PFM->cut_off=PFM->cut_off_init;
	double drad=(1.0+PFM->redshift_start)/100.0;
	for(int i=0;i<101;i++)
	  {
	    PF->rad[i]=i*drad+1.0;
	    PF->grow[i]=Growth(PFM->omega_start,PFM->lambda_start,1.0/PF->rad[i]-1.0);
	    FF << i << " " << PF->rad[i] << " " << PF->grow[i] << "\n";
	  }
      }
    int count=0;
    make_particles(fractal_memory,fractal,count,m,false);
    FF << "size " << count << "\n";
    PFM->number_particles=count;
    PF->set_number_particles(count);
    update_rv(fractal,0,0.0,0.0);
    FF << "make parts " << count << "\n";
    double delta_z=Growth(PFM->omega_0,PFM->omega_lambda,PFM->redshift_start);
    double vfratio=1.0/(1.5*PFM->omega_start)*dGrowthdT(PFM->omega_start,PFM->lambda_start,0.0);
    double omega_fraction=1.0/(1.5*PFM->omega_start);
    PF->omega_fraction=omega_fraction;
    PFM->make_scaling();
******


******

  }
}
