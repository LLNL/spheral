#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::particles="cosmos_zel";
  //
  template <class M, class F>  void make_particles(M& mem,F& frac,int& count,const double& m,const bool& crash)
  {
    ofstream& FileFractal=mem.p_fractal->p_file->DUMPS;
    //    ofstream& FileFractal=mem.p_fractal->p_file->FileFractal;
    double rand_max=RAND_MAX;
    int length=mem.grid_length;
    int length_2=length*length;
    assert(length*length_2==mem.number_particles);
    double delta=1.0/(double)length;
    if(!crash)
      {
	frac.particle_list.resize(mem.number_particles);
	for(int n=0;n<mem.number_particles;++n)
	  {
	    double x0,y0,z0;
	    if(mem.random_initial)
	      {
		x0=Fractal::my_rand(rand_max);
		y0=Fractal::my_rand(rand_max);
		z0=Fractal::my_rand(rand_max);
	      }
	    else
	      {
		x0=((double)(n % length)+mem.off_x)*delta;
		y0=((double)(n/length % length)+mem.off_y)*delta;
		z0=((double)(n/length_2)+mem.off_z)*delta;
	      }
	    int many=split_particle(mem,frac,x0,y0,z0,count,m,1,true);
	    assert(many > 0);
	  }
      }
    else
      {
	int iimax=100;
	double crash_0=mem.redshift_start+1.0;
	double scale_crash=pow(crash_0,1.0/(double)(iimax+1));
	for(int ii=0;ii<=iimax;ii++)
	  {
	    count=0;
	    for(int n=0;n<mem.number_particles;++n)
	      {
		Particle& particle=*frac.particle_list[n];
		double x0,y0,z0;
		if(mem.random_initial)
		  {
		    x0=Fractal::my_rand(rand_max);
		    y0=Fractal::my_rand(rand_max);
		    z0=Fractal::my_rand(rand_max);
		  }
		else
		  {
		    x0=((double)(n % length)+mem.off_x)*delta;
		    y0=((double)(n/length % length)+mem.off_y)*delta;
		    z0=((double)(n/length_2)+mem.off_z)*delta;
		  }
		int split_to=1;
		if(particle.get_rad_max() < 0.0)
		  {
		    split_to=pow(-crash_0/particle.get_rad_max(),mem.crash_pow);
		    split_to=min(max(1,split_to),mem.crash_levels);
		  }
		int many=split_particle(mem,frac,x0,y0,z0,count,m,split_to,false);
		assert(many > 0);
		if(count > mem.max_particles) break;
	      }
	    FileFractal << "iimax " << ii << " " << crash_0 << " " << count << :\n";
	    if(count <= mem.max_particles) break;
	    crash_0/=scale_crash;
	  }
	vector <double> crash(mem.number_particles);
	for(int n=0;n<mem.number_particles;++n)
	  {
	    Particle& particle=*frac.particle_list[n];
	    crash[n]=particle.get_rad_max();
	    delete frac.particle_list[n];
	    frac.particle_list[n]=0;
	  }
	//
	frac.particle_list.resize(count);
	frac.set_number_particles(count);
	count=0;
	for(int n=0;n<mem.number_particles;++n)
	  {
	    double x0,y0,z0;
	    if(mem.random_initial)
	      {
		x0=Fractal::my_rand(rand_max);
		y0=Fractal::my_rand(rand_max);
		z0=Fractal::my_rand(rand_max);
	      }
	    else
	      {
		x0=((double)(n % length)+mem.off_x)*delta;
		y0=((double)(n/length % length)+mem.off_y)*delta;
		z0=((double)(n/length_2)+mem.off_z)*delta;
	      }
	    int split_to=1;
	    if(crash[n] < 0.0)
	      {
		split_to=pow(-crash_0/crash[n],mem.crash_pow);
		split_to=min(max(1,split_to),mem.crash_levels);
	      }
	    int many=split_particle(mem,frac,x0,y0,z0,count,m,split_to,true);
	    assert(many > 0);
	    assert(count < mem.max_particles);
	  }
	crash.clear();
	FileFractal << "startit " << count << " " << crash_0 << " " << mem.crash_levels << :\n";
	mem.number_particles=count;
      }
  }
  template void make_particles(Fractal_Memory& mem,Fractal& frac,int& count,const double& m,const bool& crash);
}
