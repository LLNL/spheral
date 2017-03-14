#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void test_gal(Fractal_Memory& mem,Fractal& fractal)
  {
    ofstream& FileFractal=mem.p_file->DUMPS;
    //    ofstream& FileFractal=mem.p_file->FileFractal;
    FileFractal << " Made It test_gal a " << "\n";
    //    double rand_max=(double)RAND_MAX;
    //    double rmax=1.0e-5;
    double rmax=0.3;
    double x_off=0.405;
    double y_off=0.46;
    double z_off=0.4247;
    double total_mass=1.0;
    double pottrue=0.0;
    double ftrue=0.0;
    vector <double>pos(3);
    vector <double>field(4);
    double eps=pow(2.0,fractal.get_level_max())*static_cast<double>(fractal.get_grid_length());
    double eps2=1.0/(eps*eps);
    for(int n=0;n < fractal.get_number_particles();++n)
      {
	//	FileFractal << "fractal " << n << "\n";
	Particle* p=fractal.particle_list[n];
	if(p->get_p_highest_level_group() == 0)
	  continue;
	//	double m=p->get_mass();
	p->get_pos(pos);
	p->get_field_pf(field);
	double dx=pos[0]-x_off;
	double dy=pos[1]-y_off;
	double dz=pos[2]-z_off;
	double dr2=dx*dx+dy*dy+dz*dz+1.0e-30;
	double dr=sqrt(dr2);
	double fr=(dx*field[1]+dy*field[2]+dz*field[3])/dr;
	double frx=fr*dx/dr;
	double fry=fr*dy/dr;
	double frz=fr*dz/dr;
	double ft=sqrt(pow(field[1]-frx,2)+pow(field[2]-fry,2)+pow(field[3]-frz,2));
	if(dr < rmax)
	  {
	    pottrue=total_mass/rmax*(log(dr/rmax)-1.0);
	    ftrue=-total_mass*dr/((dr2+eps2)*rmax);
	  }
	else
	  {
	    pottrue=-total_mass/dr;
	    ftrue=-total_mass/(dr*dr);
	  }
	FileFractal << " test a " << n << scientific << "\t" << dr << "\t";
	FileFractal << -field[0] << "\t" << -pottrue << "\t";
	FileFractal << -fr << "\t" << -ftrue << "\t" << ft << "\t";
	FileFractal << field[1] << "\t" << frx << "\t";
	FileFractal << field[2] << "\t" << fry << "\t";
	FileFractal << field[3] << "\t" << frz << "\t";
	FileFractal << " l" << p->get_highest_level() << "\n";
      }
    FileFractal << " Made It test_gal b " << "\n";
  }
}
