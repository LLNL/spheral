#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int fractal_force_wrapper(Fractal_Memory* PFM,Fractal* PF)
  {
    cout << "starting wrapper " << endl;
    //    fractal_force_init(PFM,PF);
    PFM->p_fractal=PF;
    ofstream& FileFractal=PFM->p_file->FileFractal;
    FileFractal << " hahaa " << PFM << " " << PF << endl;
    Fractal_Memory& FM=*PFM;
    FileFractal << " hahab " << PFM<< endl;
    Fractal& FR=*PF;
    FileFractal << " hahac " << PF << " " << PF->p_file << endl;
    FileFractal << FM.p_mess->FractalRank << " starting out " << endl;
    vector <int> Box(6);
    FR.getBox(Box);
    FileFractal << " Box wrapper ";
    Misc::vector_print(Box,FileFractal);
    //
    srand(FM.random_gen+FM.p_mess->FractalRank);
    double pi=atan(1.0)*4.0;
    if(FM.periodic)
      {
	double m=0.375*FM.omega_start/pi/(double)FM.p_mess->number_particles_total;
	FileFractal << "m= " << m << endl;
	PFM->base_mass=m;
	PF->set_base_mass(m);
	bool zel_predict=FM.crash_levels > 0 && FM.max_particles > FM.number_particles; 
	int splits_tmp=0;
	double cut_off_tmp=0.0;
	if(zel_predict)
	  {
	    splits_tmp=FM.splits;
	    FM.splits=0;
	    cut_off_tmp=FM.cut_off;
	    FM.cut_off=FM.cut_off_init;
	    double drad=(1.0+FM.redshift_start)/100.0;
	    for(int i=0;i<101;i++)
	      {
		FR.rad[i]=i*drad+1.0;
		FR.grow[i]=Growth(FM.omega_start,FM.lambda_start,1.0/FR.rad[i]-1.0);
		FileFractal << i << " " << FR.rad[i] << " " << FR.grow[i] << endl;
	      }
	  }
	int count=0;
	make_particles(FM,FR,count,m,false);
	FileFractal << "size " << count << endl;
	FM.number_particles=count;
	FR.set_number_particles(count);
	update_rv(FR,0,0.0,0.0);
	FileFractal << "make parts " << count << endl;
	double delta_z=Growth(FM.omega_0,FM.omega_lambda,FM.redshift_start);
	double vfratio=1.0/(1.5*FM.omega_start)*dGrowthdT(FM.omega_start,FM.lambda_start,0.0);
	double omega_fraction=1.0/(1.5*FM.omega_start);
	FR.omega_fraction=omega_fraction;
	FM.make_scaling();
	FR.timing(-2,0);
	FR.timing(-1,49);
	fractal_force(FR,FM);
	FR.timing(1,49);
	FR.timing(0,0);
	FR.timing_lev(0,0);
	FileFractal << "initial " << delta_z << " " << vfratio << " " << omega_fraction << endl;
	if(zel_predict)
	  {
	    FM.splits=splits_tmp;
	    FM.cut_off=cut_off_tmp;
	    srand(FM.random_gen+FM.p_mess->FractalRank);
	    int count=0;
	    make_particles(FM,FR,count,m,true);
	    FileFractal << "size " << count << endl;
	    FM.number_particles=count;
	    FR.set_number_particles(count);
	    update_rv(FR,0,0.0,0.0);
	    FileFractal << "make parts " << count << endl;
	    double delta_z=Growth(FM.omega_0,FM.omega_lambda,FM.redshift_start);
	    double vfratio=1.0/(1.5*FM.omega_start)*dGrowthdT(FM.omega_start,FM.lambda_start,0.0);
	    double omega_fraction=1.0/(1.5*FM.omega_start);
	    FR.omega_fraction=omega_fraction;
	    FM.make_scaling();
	    FR.timing(-2,0);
	    FR.timing(-1,49);
	    fractal_force(FR,FM);
	    FR.timing(1,49);
	    FR.timing(0,0);
	    FR.timing_lev(0,0);
	    FileFractal << "initial " << delta_z << " " << vfratio << " " << omega_fraction << endl;
	    //	    write_rv(-6,FR);
	  }
	if(FM.norm_what == 1)omega_fraction/=vfratio;
	update_rv(FR,1,omega_fraction,0.0);
	FM.start_up=false;
	FM.calc_density_particle=false;
	FM.calc_shear=false;
	FR.timing(-2,0);
	FR.timing(-1,49);
	fractal_force(FR,FM);
	FR.timing(1,49);
	FR.timing(0,0);
	FR.timing_lev(0,0);
	//	write_rv(-5,FR);
	double lambda=Lambda(FM.omega_0,FM.omega_lambda,FM.redshift_start);
	double dp=-FM.step_length*0.5;
	if(Fractal::integrator != "leapfrog") dp=0.0;
	double omega_fraction_v=dGrowthdT(FM.omega_start,lambda,1.0/(1.0+dp)-1.0);
	update_rv(FR,2,omega_fraction,omega_fraction*omega_fraction_v);
	//	write_rv(-4,FR);
	FM.arad=1.0;
	FM.time=Age_of_the_universe(FM.omega_start,FM.lambda_start,0.0);
      }
    else
      {
	int count=0;
	FileFractal << " making particles a " << endl;
	double m=FM.total_mass/FM.p_mess->number_particles_total;      
	FileFractal << " making particles b " << m << endl;
	make_particles(FM,FR,count,m,false);
	FileFractal << "size " << count << endl;
	FM.number_particles=count;
	FR.set_number_particles(count);
	FileFractal << "make parts " << count << endl;
	FM.start_up=false;
	FM.calc_density_particle=false;
	FM.calc_shear=false;
	FR.timing(-2,0);
	FR.timing(-1,49);
	fractal_force(FR,FM);
	FR.timing(1,49);
	FR.timing(0,0);
	FR.timing_lev(0,0);
	velocities(FM,FR);
	//	write_rv(-3,FR);
      }
    //
    if(FM.number_steps_total == 0)
      {
	int iphase=3;
	int jfield=4;
	if(FM.calc_density_particle) jfield=5;
	//	if(FM.calc_shear) jfield=12;
	fix_memory(FR,iphase,jfield);
	FR.timing(-2,0);
	FR.timing(-1,49);
	fractal_force(FR,FM);
	FR.timing(1,49);
	FR.timing(0,0);
	FR.timing_lev(0,0);
	return 0;
      }
    for(int step=0;step < FM.number_steps_total; ++step)
      {
	if(FM.steps == -1) energy_simple(FM,FR);
	FM.steps=step;
	FileFractal << "step = " << step << " " << FM.arad << endl;
	if(FM.periodic && step % FM.number_steps_out == 0)
	  {
	    FM.calc_density_particle=true;
	    FM.do_var=true;
	  }
	else
	  {
	    FM.calc_density_particle=false;
	    FM.do_var=false;
	  }
	int iphase=6;
	int jfield=4;
	if(FM.calc_density_particle) jfield=5;
	fix_memory(FR,iphase,jfield);
	step_simple(FM,FR);
	energy_simple(FM,FR);
	if(step % FM.number_steps_out == 0)
	  write_rv(step,FR);
      }
    return 0;
  }
}
