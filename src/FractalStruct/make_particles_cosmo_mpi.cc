#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::particles="cosmos_zel";
  //
  template <class M, class F>  void make_particles(M& mem,F& frac,int& count,const double& m,const bool& crashcrash)
  {
    ofstream& FileFractal=mem.p_fractal->p_file->DUMPS;
    //    ofstream& FileFractal=mem.p_fractal->p_file->FileFractal;
    vector <int> Box(6);
    frac.getBox(Box);
    vector <int> BoxLength(3);
    BoxLength[0]=Box[1]-Box[0]+1;
    BoxLength[1]=Box[3]-Box[2]+1;
    BoxLength[2]=Box[5]-Box[4]+1;
    int total=BoxLength[0]*BoxLength[1]*BoxLength[2];
    mem.number_particles=total;
    int length=mem.grid_length;
    double delta=1.0/(double)length;
    double x0,y0,z0;
    if(!crashcrash)
      {
	FileFractal << " Boxfrac " << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << " " << Box[4] << " " << Box[5] << "\n";
	frac.particle_list.resize(mem.number_particles);
	for(int nx=Box[0];nx<=Box[1];nx++)
	  {
	    x0=(static_cast<double>(nx)+mem.off_x)*delta;
	    for(int ny=Box[2];ny<=Box[3];ny++)
	      {
		y0=(static_cast<double>(ny)+mem.off_y)*delta;
		for(int nz=Box[4];nz<=Box[5];nz++)
		  {
		    z0=(static_cast<double>(nz)+mem.off_z)*delta;
		    //		    FileFractal << " makeb " << nx << " " << ny << " " << nz << " " << count << "\n";
		    int many=split_particle(mem,frac,x0,y0,z0,count,m,1,true);
		    //		    FileFractal << " makec " << nx << " " << ny << " " << nz << " " << many << " " << count << "\n";
		    assert(many > 0);
		  }
	      }
	  }
      }
    else
      {
	int iimax=400;
	double crash_0=mem.redshift_start+1.0;
	double scale_crash=pow(crash_0,2.0/(double)(iimax+1));
	bool breakA=false;
	for(int ii=0;ii<=iimax;ii++)
	  {
	    int n=0;
	    count=0;
	    for(int nx=Box[0];nx<=Box[1];nx++)
	      {
		x0=(static_cast<double>(nx)+mem.off_x)*delta;
		for(int ny=Box[2];ny<=Box[3];ny++)
		  {
		    y0=(static_cast<double>(ny)+mem.off_y)*delta;
		    for(int nz=Box[4];nz<=Box[5];nz++)
		      {
			z0=(static_cast<double>(nz)+mem.off_z)*delta;
			Particle& particle=*frac.particle_list[n];
			int split_to=1;
			if(particle.get_real_particle() && particle.get_rad_max() < 0.0)
			  {
			    split_to=pow(-crash_0/particle.get_rad_max(),mem.crash_pow);
			    split_to=min(max(1,split_to),mem.crash_levels);
			  }
			n++;
			int many=split_particle(mem,frac,x0,y0,z0,count,m,split_to,false);
			//			FileFractal << " RAD " << crash_0 << " " << mem.max_particles << " " << n-1 << " " << count << " " << particle.get_rad_max() << " " << split_to << "\n";
			assert(many > 0);
			breakA = count > mem.max_particles;
			if(breakA)
			  break;
		      }
		    if(breakA)
		      break;
		  }
		if(breakA)
		  break;
	      }
	    FileFractal << "iimax " << ii << " " << crash_0 << " " << count << " " << n << "\n";
	    if(!breakA)
	      break;
	    crash_0/=scale_crash;
	  }
	FileFractal << "iimax " << crash_0<< "\n";;
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
	double x0,y0,z0;
	count=0;
	int n=0;
	for(int nx=Box[0];nx<=Box[1];nx++)
	  {
	    x0=(static_cast<double>(nx)+mem.off_x)*delta;
	    for(int ny=Box[2];ny<=Box[3];ny++)
	      {
		y0=(static_cast<double>(ny)+mem.off_y)*delta;
		for(int nz=Box[4];nz<=Box[5];nz++)
		  {
		    z0=(static_cast<double>(nz)+mem.off_z)*delta;
		    int split_to=1;
		    if(crash[n] < 0.0)
		      {
			split_to=pow(-crash_0/crash[n],mem.crash_pow);
			split_to=min(max(1,split_to),mem.crash_levels);
		      }
		    int many=split_particle(mem,frac,x0,y0,z0,count,m,split_to,true);
		    assert(many > 0);
		    if(count >= mem.max_particles)
		      {
			FileFractal << " OVER " << count << " " << mem.max_particles;
			FileFractal << " " << Box[0] << " " << nx << " " << Box[1];
			FileFractal << " " << Box[2] << " " << ny << " " << Box[3];
			FileFractal << " " << Box[4] << " " << nz << " " << Box[5] << "\n";
		      }
		    //		    assert(count < mem.max_particles);
		    n++;
		  }
	      }
	  }
	crash.clear();
	FileFractal << "startit " << count << " " << mem.max_particles << " " << crash_0 << " " << mem.crash_levels << "\n";
	mem.number_particles=count;
      }
  }
}
namespace FractalSpace
{
  template void make_particles(Fractal_Memory& mem,Fractal& frac,int& count,const double& m,const bool& crashcrash);
}
