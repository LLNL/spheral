#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int fractal_force_wrapper(Fractal_Memory* p_fractal_memory,Fractal* p_fractal)
  {
    cout << "starting wrapper " << endl;
    //    fractal_force_init(p_fractal_memory,p_fractal);
    p_fractal_memory->p_fractal=p_fractal;
    ofstream& FileFractal=p_fractal_memory->p_file->FileFractal;
    FileFractal << " hahaa " << p_fractal_memory << " " << p_fractal << endl;
    Fractal_Memory& fractal_memory=*p_fractal_memory;
    FileFractal << " hahab " << p_fractal_memory<< endl;
    Fractal& fractal=*p_fractal;
    FileFractal << " hahac " << p_fractal << " " << p_fractal->p_file << endl;
    FileFractal << fractal_memory.p_mess->FractalRank << " starting out " << endl;
    vector <int> Box(6);
    fractal.getBox(Box);
    FileFractal << " Box wrapper ";
    Misc::vector_print(Box,FileFractal);
    //
    srand(fractal_memory.random_gen+fractal_memory.p_mess->FractalRank);
    double pi=atan(1.0)*4.0;
    if(fractal_memory.periodic)
      {
	double m=0.375*fractal_memory.omega_start/pi/(double)fractal_memory.p_mess->number_particles_total;
	FileFractal << "m= " << m << endl;
	bool zel_predict=fractal_memory.crash_levels > 0 && fractal_memory.max_particles > fractal_memory.number_particles; 
	int splits_tmp=0;
	double cut_off_tmp=0.0;
	if(zel_predict)
	  {
	    splits_tmp=fractal_memory.splits;
	    fractal_memory.splits=0;
	    cut_off_tmp=fractal_memory.cut_off;
	    fractal_memory.cut_off=fractal_memory.cut_off_init;
	    double drad=(1.0+fractal_memory.redshift_start)/100.0;
	    for(int i=0;i<101;i++)
	      {
		fractal.rad[i]=i*drad+1.0;
		fractal.grow[i]=Growth(fractal_memory.omega_start,fractal_memory.lambda_start,1.0/fractal.rad[i]-1.0);
		FileFractal << i << " " << fractal.rad[i] << " " << fractal.grow[i] << endl;
	      }
	  }
	int count=0;
	make_particles(fractal_memory,fractal,count,m,false);
	FileFractal << "size " << count << endl;
	fractal_memory.number_particles=count;
	fractal.set_number_particles(count);
	update_rv(fractal,0,0.0,0.0);
	FileFractal << "make parts " << count << endl;
	double delta_z=Growth(fractal_memory.omega_0,fractal_memory.omega_lambda,fractal_memory.redshift_start);
	double vfratio=1.0/(1.5*fractal_memory.omega_start)*dGrowthdT(fractal_memory.omega_start,fractal_memory.lambda_start,0.0);
	double omega_fraction=1.0/(1.5*fractal_memory.omega_start);
	fractal.omega_fraction=omega_fraction;
	fractal_memory.make_scaling();
	fractal.timing(-2,0);
	fractal.timing(-1,49);
	fractal_force(fractal,fractal_memory);
	fractal.timing(1,49);
	fractal.timing(0,0);
	fractal.timing_lev(0,0);
	FileFractal << "initial " << delta_z << " " << vfratio << " " << omega_fraction << endl;

	if(zel_predict)
	  {
	    fractal_memory.splits=splits_tmp;
	    fractal_memory.cut_off=cut_off_tmp;
	    int count=0;
	    make_particles(fractal_memory,fractal,count,m,true);
	    FileFractal << "size " << count << endl;
	    fractal_memory.number_particles=count;
	    fractal.set_number_particles(count);
	    update_rv(fractal,0,0.0,0.0);
	    FileFractal << "make parts " << count << endl;
	    double delta_z=Growth(fractal_memory.omega_0,fractal_memory.omega_lambda,fractal_memory.redshift_start);
	    double vfratio=1.0/(1.5*fractal_memory.omega_start)*dGrowthdT(fractal_memory.omega_start,fractal_memory.lambda_start,0.0);
	    double omega_fraction=1.0/(1.5*fractal_memory.omega_start);
	    fractal.omega_fraction=omega_fraction;
	    fractal_memory.make_scaling();
	    fractal.timing(-2,0);
	    fractal.timing(-1,49);
	    fractal_force(fractal,fractal_memory);
	    fractal.timing(1,49);
	    fractal.timing(0,0);
	    fractal.timing_lev(0,0);
	    FileFractal << "initial " << delta_z << " " << vfratio << " " << omega_fraction << endl;
	  }
	if(fractal_memory.norm_what == 1)omega_fraction/=vfratio;
	update_rv(fractal,1,omega_fraction,0.0);
	fractal_memory.start_up=false;
	fractal_memory.calc_density_particle=false;
	fractal_memory.calc_shear=false;
	fractal.timing(-2,0);
	fractal.timing(-1,49);
	fractal_force(fractal,fractal_memory);
	fractal.timing(1,49);
	fractal.timing(0,0);
	fractal.timing_lev(0,0);
	double lambda=Lambda(fractal_memory.omega_0,fractal_memory.omega_lambda,fractal_memory.redshift_start);
	double dp=-fractal_memory.step_length*0.5;
	if(Fractal::integrator != "leapfrog") dp=0.0;
	double omega_fraction_v=dGrowthdT(fractal_memory.omega_start,lambda,1.0/(1.0+dp)-1.0);
	update_rv(fractal,2,omega_fraction,omega_fraction*omega_fraction_v);
	fractal_memory.arad=1.0;
	fractal_memory.time=Age_of_the_universe(fractal_memory.omega_start,fractal_memory.lambda_start,0.0);
      }
    else
      {
	int count=0;
	FileFractal << " making particles a " << endl;
	double m=fractal_memory.total_mass/fractal_memory.p_mess->number_particles_total;      
	FileFractal << " making particles b " << m << endl;
	make_particles(fractal_memory,fractal,count,m,false);
	FileFractal << "size " << count << endl;
	fractal_memory.number_particles=count;
	fractal.set_number_particles(count);
	FileFractal << "make parts " << count << endl;
	fractal_memory.start_up=false;
	fractal_memory.calc_density_particle=false;
	fractal_memory.calc_shear=false;
	fractal.timing(-2,0);
	fractal.timing(-1,49);
	fractal_force(fractal,fractal_memory);
	fractal.timing(1,49);
	fractal.timing(0,0);
	fractal.timing_lev(0,0);
	velocities(fractal_memory,fractal);
      }
    //
    if(fractal_memory.number_steps_total == 0)
      {
	int iphase=3;
	int jfield=4;
	if(fractal_memory.calc_density_particle) jfield=5;
	//	if(fractal_memory.calc_shear) jfield=12;
	fix_memory(fractal,iphase,jfield);
	fractal.timing(-2,0);
	fractal.timing(-1,49);
	fractal_force(fractal,fractal_memory);
	fractal.timing(1,49);
	fractal.timing(0,0);
	fractal.timing_lev(0,0);
	return 1;
      }
    for(int step=0;step < fractal_memory.number_steps_total; ++step)
      {
	if(fractal_memory.steps == -1) energy_simple(fractal_memory,fractal);
	fractal_memory.steps=step;
	FileFractal << "step = " << step << " " << fractal_memory.arad << endl;
	if(fractal_memory.periodic && step % fractal_memory.number_steps_out == 0)
	  {
	    fractal_memory.calc_density_particle=true;
	    fractal_memory.do_var=true;
	  }
	else
	  {
	    fractal_memory.calc_density_particle=false;
	    fractal_memory.do_var=false;
	  }
	int iphase=6;
	int jfield=4;
	if(fractal_memory.calc_density_particle) jfield=5;
	fix_memory(fractal,iphase,jfield);
	step_simple(fractal_memory,fractal);
	energy_simple(fractal_memory,fractal);
	if(step % fractal_memory.number_steps_out == 0)
	  write_rv(step,fractal);
      }
    return 1;
  }
}
