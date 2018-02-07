#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_pseudo_particles(Fractal_Memory& mem,Fractal& frac)
  {
    //    frac.p_file->FileFractal << " entered into add pseudo particles " << "\n";
    if(!mem.periodic)
      return;
    //    ofstream& FF=mem.p_file->FileFractal;
    int length=mem.grid_length;
    double Rdelta=1.0/static_cast<double>(length);
    double Rlow=-2.0*Rdelta;
    double Rhigh=1.0+Rdelta;
    vector <double>pos(3);
    // vector <double> boxinner(6);
    // vector <double> boxouter(6);
    // boxinner[0]=0.0;
    // boxinner[1]=1.0;
    // boxinner[2]=0.0;
    // boxinner[3]=1.0;
    // boxinner[4]=0.0;
    // boxinner[5]=1.0;
    // boxouter[0]=Rlow;
    // boxouter[1]=Rhigh;
    // boxouter[2]=Rlow;
    // boxouter[3]=Rhigh;
    // boxouter[4]=Rlow;
    // boxouter[5]=Rhigh;
    vector<double>boxinner{0.0,1.0,0.0,1.0,0.0,1.0};
    vector<double>boxouter{Rlow,Rhigh,Rlow,Rhigh,Rlow,Rhigh};
    clean_deque(frac.pseudo_particle_list);
    vector <double> posp(3);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	//	FF << " add " << particle << " " ;
	//	P->dump(mem.p_file->FileFractal);
	P->get_pos(pos);
	for(int nz=-1;nz<=1;nz++)
	  {
	    posp[2]=pos[2]+nz;
	    for(int ny=-1;ny<=1;ny++)
	      {
		posp[1]=pos[1]+ny;
		for(int nx=-1;nx<=1;nx++)
		  {
		    posp[0]=pos[0]+nx;
		    if(nz==0 && ny==0 && nx==0)
		      continue;
		    if(vector_in_box(posp,boxinner))
		       continue;
		    if(!vector_in_box(posp,boxouter))
		       continue;
		    double m=P->get_mass();
		    Particle* Pb=new Particle;
		    Pb->set_posm(posp,m);
		    Pb->set_real_particle(false);
		    frac.particle_list.push_back(Pb);
		    frac.pseudo_particle_list.push_back(Pb);
		    //		    Pb->dump(mem.p_file->FileFractal);
		  }
	      }
	  }
      }
    int parts=frac.particle_list.size();
    frac.set_number_particles(parts);
  }
}
