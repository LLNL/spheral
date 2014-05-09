#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::energy_method("leapfrog");
  //
  template <class M>  void energy_simple(M& mem,Fractal& fractal)
  {
    ofstream& FileEnergy=fractal.p_file->FileEnergy;
    ofstream& FileMom=fractal.p_file->DUMPS;
    //    ofstream& FileMom=fractal.p_file->FileMom;
    ofstream& FileP=fractal.p_file->DUMPS;
    //    ofstream& FileP=fractal.p_file->FileParticle;
    assert(Fractal::integrator=="leapfrog");
    double sum_m=0.0;
    double sum_x=0.0;
    double sum_y=0.0;
    double sum_z=0.0;
    double sum_vx=0.0;
    double sum_vy=0.0;
    double sum_vz=0.0;
    double sum_ekx=0.0;
    double sum_eky=0.0;
    double sum_ekz=0.0;
    double sum_p=0.0;
    vector <double>pos(3,0.0);
    vector <double>vel(3,0.0);
    FileP << " number of particles " << fractal.get_number_particles() << "\n";
    for(int n=0;n < fractal.get_number_particles();++n)
      {
	Particle* p=fractal.particle_list[n];
	if(p->get_p_highest_level_group() == 0)
	  continue;
	double m=p->get_mass();
	p->get_phase(pos,vel);
	sum_m+=m;
	sum_x+=m*pos[0];
	sum_y+=m*pos[1];
	sum_z+=m*pos[2];
	sum_vx+=m*vel[0];
	sum_vy+=m*vel[1];
	sum_vz+=m*vel[2];
	sum_ekx+=m*vel[0]*vel[0];
	sum_eky+=m*vel[1]*vel[1];
	sum_ekz+=m*vel[2]*vel[2];
	sum_p+=m*p->get_potential();
      }
    sum_x/=sum_m;
    sum_y/=sum_m;
    sum_z/=sum_m;
    double sum_px=0.0;
    double sum_py=0.0;
    double sum_pz=0.0;
    for(int n=0;n < fractal.get_number_particles();++n)
      {
	Particle* p=fractal.particle_list[n];
	if(p->get_p_highest_level_group() == 0)
	  continue;
	double m=p->get_mass();
	p->get_phase(pos,vel);
	sum_px=m*(pos[1]*vel[2]-pos[2]*vel[1]);
	sum_py=m*(pos[2]*vel[0]-pos[0]*vel[2]);
	sum_pz=m*(pos[0]*vel[1]-pos[1]*vel[0]);
      }
    FileP << " sum " << sum_vx << " " << sum_vy << " " << sum_vz << " " << sum_ekx << " " << sum_eky << " " << sum_ekz << " " << sum_p << "\n";
    mem.potential_energy=0.5*sum_p;
    mem.kinetic_energy=0.5*(sum_ekx+sum_eky+sum_ekz);
    vector <double>sums(9);
    if(mem.periodic)
      {
	double parad=pow(mem.arad,mem.pexp);
	if(mem.steps == 0)
	  {
	    double arad_down=pow(parad-0.5*mem.step_length,1.0/mem.pexp);
	    mem.potential_energy_old=mem.potential_energy*arad_down*arad_down;
	    mem.udda=-mem.potential_energy_old*arad_down/3.0;
	  }
	double parad_down=parad-0.5*mem.step_length;
	double parad_up=parad+0.5*mem.step_length;
	double arad_down=pow(parad_down,1.0/mem.pexp);
	double arad_up=pow(parad_up,1.0/mem.pexp);
	double da=arad_up-arad_down;
	double e_p=mem.potential_energy*mem.arad;
	double e_k=(mem.kinetic_energy_old*pow(arad_down,4)+
		    mem.kinetic_energy*pow(arad_up,4))/2.0;
	mem.udda-=0.5*da*(mem.potential_energy+mem.potential_energy_old);
	double total_energy=e_p+e_k+mem.udda;
	if(mem.steps >= 0)
	  {
	    FileEnergy << scientific << mem.time << "\t " << mem.arad << "\t " << mem.steps << "\t " << 
	      total_energy << "\t " << -e_p << "\t " << e_k << "\t " << mem.udda;
	    if(mem.MPIrun)
	      {
		sums[0]=e_p;
		sums[1]=e_k;
		sums[2]=mem.udda;
		sums[3]=sum_vx;
		sums[4]=sum_vy;
		sums[5]=sum_vz;
		sums[6]=sum_px;
		sums[7]=sum_py;
		sums[8]=sum_pz;
		mem.p_mess->Find_Sum_DOUBLE(sums,9);
		FileEnergy << "\t" << sums[0]+sums[1]+sums[2] << "\t" << -sums[0] << "\t" << sums[1] << "\t" << sums[2];
	      }
	    FileEnergy << "\t" << M::omega(mem.arad,mem.omega_start,mem.lambda_start);
	    FileEnergy << "\t" << M::lambda(mem.arad,mem.omega_start,mem.lambda_start);
	    FileEnergy << "\t" << M::hubble(mem.arad,mem.omega_start,mem.lambda_start) << "\n";
	    FileMom << scientific << mem.time << "\t " << mem.arad << "\t " << mem.steps << "\t " << 
	      sum_vx << "\t " << sum_vy << "\t " << sum_vz << "\t " <<
	      sum_px << "\t " << sum_py << "\t " << sum_pz;
	    if(mem.MPIrun)
	      FileMom << "\t" << sums[3] << "\t" << sums[4] << "\t" << sums[5] << "\t" << sums[6] << "\t" << sums[7] << "\t" << sums[8];
	    FileMom << "\n";
	  }
      }
    else
      {
	double total_energy=0.5*(mem.kinetic_energy_old+mem.kinetic_energy)+
	  mem.potential_energy;
	if(mem.steps >= 0)
	  {
	    FileEnergy << scientific << mem.time <<  "\t " << mem.steps << "\t " << 
	      total_energy << "\t " << -mem.potential_energy << "\t " << 0.5*(mem.kinetic_energy_old+mem.kinetic_energy) << "\t " << "\n";
	    if(mem.MPIrun)
	      {
		sums[0]=total_energy;
		sums[1]=0.5*(mem.kinetic_energy_old+mem.kinetic_energy);
		sums[2]=mem.potential_energy;
		sums[3]=sum_vx;
		sums[4]=sum_vy;
		sums[5]=sum_vz;
		sums[6]=sum_px;
		sums[7]=sum_py;
		sums[8]=sum_pz;
		mem.p_mess->Find_Sum_DOUBLE(sums,9);
		FileEnergy << "\t" << sums[0] << "\t" << sums[1] << "\t" << -sums[2] << "\t";
	      }
	    FileMom << scientific << mem.time << "\t " << mem.steps << "\t " << 
	      sum_vx << "\t " << sum_vy << "\t " << sum_vz << "\t " <<
	      sum_px << "\t " << sum_py << "\t " << sum_pz;
	    if(mem.MPIrun)
	      FileMom << "\t" << sums[3] << "\t" << sums[4] << "\t" << sums[5] << "\t" << sums[6] << "\t" << sums[7] << "\t" << sums[8];
	    FileMom << "\n";
	  }
      }
    mem.kinetic_energy_old=mem.kinetic_energy;
    mem.potential_energy_old=mem.potential_energy;
  }
}
namespace FractalSpace
{
  template void energy_simple(Fractal_Memory& mem, Fractal& fractal);
}
