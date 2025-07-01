#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  typedef ptrdiff_t pint;
  void initial_forces_sharp(Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFractal=mem.p_fractal->p_file->DUMPS;
//     int FractalRank=mem.p_mess->FractalRank;
    FileFractal << "enter initial_forces " << "\n";
    int seed=mem.random_gen+mem.p_mess->FractalRank;
    srand(seed);
    //    std::default_random_engine generator(seed);
    //    std::normal_distribution<double> distributionG(1.0,0.0);
    int highest_level_used=-1;
    for(int level=0;level <= frac.get_level_max();level++)
      {
	if(mem.all_groups[level].size() > 0)
	  highest_level_used=level;
      }
    vector <int>used(1);
    used[0]=highest_level_used;
    mem.p_mess->Find_Max_INT(used,1);
    highest_level_used=used[0];
    int length=frac.get_grid_length();
    double length_5=pow(static_cast<double>(length),-5);
    assert(length >0);
    int nyq=length/2;
    int length_c=nyq+1;
    int real_length=nyq*Misc::pow(2,highest_level_used)+1;
    vector <double> variance_rho(real_length,0.0);
    vector <double> variance_pot(real_length,0.0);
    vector <double> variance_force(real_length,0.0);
    vector <double> variance_force_s(real_length,0.0);
    double rand_max=static_cast<double>(RAND_MAX);
    double pi=4.0*atan(1.0);
    double twopi=2.0*pi;
    double fourpi=4.0*pi;
    assert(mem.highest_level_init >=0 && mem.highest_level_init <= frac.get_level_max());
    int highest_level_fft=min(mem.highest_level_init,highest_level_used);
    double var_initial=mem.sigma_initial*mem.sigma_initial;
    double dead_parrot=0.0;
    FileFractal << "mem.scaling= " << mem.scaling << " " << var_initial << "\n";
    vector <double> green_2(length+1);
    for (int k=0;k <length;++k)
      {
	double aa=pi*(double)(k)/(double)(length);
	green_2[k]=pow(2.0*sin(aa),2);
      }
    for(int lev=0;lev <= highest_level_fft;++lev)
      {
	FileFractal << "make power " << lev << " " << highest_level_fft << "\n";
	mem.p_mess->create_potC();
	double power=pow(8.0,lev);
	double scale=pow(2.0,lev);
	double force_const=fourpi/scale/scale*length_5;
	int cut_lev=-1;
	if(lev > 0) cut_lev=nyq/2;
	double step_wave=pow(2.0,lev);
	//	int begin_x=mem.p_mess->start_x;
	//	int end_x=begin_x+mem.p_mess->length_x;
	for(pint kx=mem.p_mess->start_x;kx < mem.p_mess->start_x+mem.p_mess->length_x; kx++)
	  {
	    int ka=min(kx,length-kx);
	    bool nyq_x= kx==0 || kx== nyq;
	    bool shorty_x= ka <= cut_lev;
	    for(int ky=0;ky < length ; ky++)
	      {
		int kb=min(ky,length-ky);
		bool nyq_y= ky==0 || ky== nyq;
		bool shorty_y= kb <= cut_lev;
		for(int kz=0;kz <= nyq ; kz++)
		  {
		    int kc=kz;
		    bool nyq_z= kz==0 || kc== nyq;
		    bool nyquist=nyq_x && nyq_y && nyq_z;
		    bool shorty_z= kz <= cut_lev;
		    double k=step_wave*sqrt((double)(ka*ka+kb*kb+kc*kc));
		    bool shorty=(shorty_x && shorty_y && shorty_z) || k < 0.001;
		    double amplitude_raw=0.0;
		    if(k > 0.001 && !shorty)
		      amplitude_raw=sqrt(power*cosmos_power(k/mem.scaling,mem));
		    double angle=0.0;
		    if(!nyquist)
		      angle=twopi*Fractal::my_rand(rand_max);
		    double norwegian_blue=Fractal::my_rand_not_zero(rand_max);
		    double amplitude_random=amplitude_raw*sqrt(-2.0*log(norwegian_blue));
		    double pot_k=amplitude_random*step_wave;
		    int holy_grail=mem.fftw_where(kx,ky,kz,length,length_c);
		    //
		    mem.p_mess->potC[holy_grail][0]=pot_k*cos(angle);
		    mem.p_mess->potC[holy_grail][1]=pot_k*sin(angle);
		  }
	      }
	  }
	FileFractal << "calling power_spectrum from initial_forces " << length << "\n";
	FileFractal << "sizes a " << variance_rho.size() << " " << variance_pot.size() << " " << variance_force.size() << " " << variance_force_s.size() << "\n";
	power_spectrum(mem.p_mess->potC,length,variance_rho,variance_pot,variance_force,variance_force_s,lev,frac.get_density_0(),true,mem);
	//	power_spectrum(mem.p_mess->potC,length,variance_rho,variance_pot,variance_force,variance_force_s,lev,frac.get_density_0(),true,mem);
	//
	FileFractal << "back from power " << lev << "\n";
	if(lev == 0)
	  {
	    double a=mem.norm_scale*(double)length;
	    int n1=(int)a;
	    a-=(double)n1;
	    double var_obs_0=-1.0;
	    if(mem.norm_what == 0)
	      {
		double var_obs=(1.0-a)*variance_rho[n1]+a*variance_rho[n1+1];
		dead_parrot=sqrt(var_initial/var_obs);
		FileFractal << "dead parrot " << dead_parrot << "\n";
	      }
	    else if(mem.norm_what == 1 || mem.norm_what==2)
	      {
		var_obs_0=((1.0-a)*variance_force_s[n1]+a*variance_force_s[n1+1]);
		assert(0);
	      }
	    else if(mem.norm_what == 3)
	      {
		var_obs_0=((1.0-a)*variance_pot[n1]+a*variance_pot[n1+1]);
		assert(0);
	      }
	    else
	      assert(0);
	  }
	//
	for(int i=0;i < length/2;i++)
	  {
	    if(mem.p_mess->IAmAnFFTNode)
	      FileFractal << "real variance a " << i << " " << (double)i/(double)length << " " << variance_rho[i] << "\n";
	  }
	for(pint kx=mem.p_mess->start_x;kx < mem.p_mess->start_x+mem.p_mess->length_x ; kx++)
	  {
	    FileFractal << " after var a " << mem.p_mess->FractalRank << " " << kx << "\n";
	    int ka=min(kx,length-kx);
	    for(int ky=0;ky < length ; ky++)
	      {
		int kb=min(ky,length-ky);
		for(int kz=0;kz <= nyq ; kz++)
		  {
		    int kc=kz;
		    double g2=force_const/(green_2[ka]+green_2[kb]+green_2[kc]+1.0e-30);
		    int brian=mem.fftw_where(kx,ky,kz,length,length_c);
		    mem.p_mess->potC[brian][0]*=g2;
		    mem.p_mess->potC[brian][1]*=g2;
		  }
	      }
	    FileFractal << " after var b " << mem.p_mess->FractalRank << " " << kx << "\n";
	  }
	mem.p_mess->create_potR();
	mem.p_mess->fftw_complex_to_real();
	mem.p_mess->free_potC();
	mem.p_mess->create_potRS();
	FileFractal << " POTRS " << mem.p_mess->total_memory << " " << length << " " << mem.p_mess->length_x << "\n";
	std::copy(mem.p_mess->potR,mem.p_mess->potR+2*mem.p_mess->total_memory,mem.p_mess->potRS);
	mem.p_mess->free_potR();
	Full_Stop(mem,39);
	if(!lev==0)
	  for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
	      group_itr!=mem.all_groups[lev].end();group_itr++)
	    potential_start(**group_itr);

	slices_to_potf(mem,frac,lev);
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
	    group_itr!=mem.all_groups[lev].end();group_itr++)
	  force_at_point(**group_itr,frac);
	FileFractal << " Finished LEVEL " << lev << "\n";
      }
    if(highest_level_fft < highest_level_used)
      {
	for(int lev=highest_level_fft+1;lev <= highest_level_used;++lev)
	  {
	    FileFractal << " Starting A LEVEL " << lev << "\n";
	    for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
		group_itr!=mem.all_groups[lev].end();group_itr++)
	      {
		Group& group=**group_itr;
		potential_start(group);
		force_at_point(group,frac);
	      } 
	    FileFractal << " Finishing A LEVEL " << lev << "\n";
	  }
      }
    for(int level=0;level <= frac.get_level_max();level++)
      {
	FileFractal << " Starting B LEVEL " << level << "\n";
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	    group_itr!=mem.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    group.scale_pot_forces(dead_parrot);
	    double varx,vary,varz;
	    group.get_force_variance(varx,vary,varz);
	    //	    FileFractal << "point variances " << group.get_level() << " " << group.list_points.size() << " " << varx << " " << vary << " " << varz << "\n";
	    //	    FileFractal << "sharpiea " << frac.get_level_max() << "\n";
	    force_at_particle_sharp(group,frac); 		
	    //	    FileFractal << "sharpieb " << frac.get_level_max() << "\n";
	  }
	FileFractal << " Finishing B LEVEL " << level << "\n";
      }
    FileFractal << "calling fractal " << &frac << "\n";
    sum_pot_forces(frac);
    //  assert(0);
    FileFractal << "whatt0 " << "\n";
    for(int lev=0;lev<=highest_level_used;lev++)
      {
	double sum0=1.0e-10;
	vector <double>sum1(3,1.0e-10);
	vector <double>sum2(3,1.0e-10);
	for(int i=0;i<frac.get_number_particles();i++)
	  {
	    Particle& p=*frac.particle_list[i];
	    if(!p.get_real_particle())
	      continue;
	    if((p.get_p_highest_level_group())->get_level() == lev)
	      {
		vector <double>force(3);
		p.get_force(force);
		sum0+=1.0;
		sum1[0]+=force[0];
		sum1[1]+=force[1];
		sum1[2]+=force[2];
		sum2[0]+=force[0]*force[0];
		sum2[1]+=force[1]*force[1];
		sum2[2]+=force[2]*force[2];
	      }
	  }
	sum1[0]/=sum0;
	sum1[1]/=sum0;
	sum1[2]/=sum0;
	sum2[0]/=sum0;
	sum2[1]/=sum0;
	sum2[2]/=sum0;
	sum2[0]=sqrt(sum2[0]-sum1[0]*sum1[0]);
	sum2[1]=sqrt(sum2[1]-sum1[1]*sum1[1]);
	sum2[2]=sqrt(sum2[2]-sum1[2]*sum1[2]);
	FileFractal << " level forces " << lev  << " " << sum1[0] << " " << sum2[0] << " " << sum1[1] << " " << sum2[1] << " " << sum1[2] << " " << sum2[2] << "\n";
      } 
    sum_pot_forces(frac);
    //
    //    mem.p_mess->free_potRC();
  }
}
