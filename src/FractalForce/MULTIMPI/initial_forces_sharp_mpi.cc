#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void initial_forces_sharp(Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFractal=mem.p_fractal->p_file->FileFractal;
    FileFractal << "enter initial_forces " << endl;
    int highest_level_used=-1;
    for(int level=0;level <= frac.get_level_max();level++)
      {
	if(mem.all_groups[level].size() > 0)
	  highest_level_used=level;
      }
    int* used= new int[1];
    used[0]=highest_level_used;
    mem.p_mess->Find_Max_INT(used,1);
    highest_level_used=used[0];
    delete [] used;
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
    mem.p_mess->create_potRC();
    double pi=4.0*atan(1.0);
    double twopi=2.0*pi;
    double fourpi=4.0*pi;
    assert(mem.highest_level_init >=0 && mem.highest_level_init <= frac.get_level_max());
    int highest_level_fft=min(mem.highest_level_init,highest_level_used);
    double var_initial=mem.sigma_initial*mem.sigma_initial;
    double dead_parrot=0.0;
    FileFractal << "mem.scaling= " << mem.scaling << " " << var_initial << endl;
    vector <double> green_2(length+1);
    for (int k=0;k <length;++k)
      {
	double aa=pi*(double)(k)/(double)(length);
	green_2[k]=pow(2.0*sin(aa),2);
      }
    for(int lev=0;lev <= highest_level_fft;++lev)
      {
	FileFractal << "make power " << lev << " " << highest_level_fft << endl;
	double boost_power=pow(8.0,lev);
	double boost_scale=pow(2.0,lev);
	double force_const=fourpi/boost_scale/boost_scale*length_5;
	int cut_lev=-1;
	if(lev > 0) cut_lev=nyq/2;
	double step_wave=pow(2.0,lev);
	//	int begin_x=mem.p_mess->start_x;
	//	int end_x=begin_x+mem.p_mess->length_x;
	for(ptrdiff_t kx=mem.p_mess->start_x;kx < mem.p_mess->start_x+mem.p_mess->length_x; kx++)
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
		      amplitude_raw=sqrt(boost_power*cosmos_power(k/mem.scaling,mem));
		    double angle=0.0;
		    if(nyquist)
		      {
			angle=0.0;
			amplitude_raw*=0.5;
		      }
		    else
		      {
			angle=twopi*Fractal::my_rand(rand_max);
		      }
		    //
		    //		    if(kz != 1 || ky != 0 || kx != 0) 
		    //		      amplitude_raw=1.0e-30;
		    //		    else
		    //		      angle=pi*0.5;
		    //
		    double norwegian_blue=Fractal::my_rand_not_zero(rand_max);
		    double amplitude_random=amplitude_raw*sqrt(-2.0*log(norwegian_blue));
		    double pot_k=amplitude_random*step_wave;
		    int holy_grail=mem.fftw_where(kx,ky,kz,length,length_c);
		    //
		    mem.p_mess->potC[holy_grail][0]=pot_k*cos(angle);
		    mem.p_mess->potC[holy_grail][1]=pot_k*sin(angle);
		    //		    FileFractal << " power " << aa << " "  << bb << " "  << cc << " "  << dd << " " ;
		    //		    FileFractal << mem.p_mess->potC[holy_grail][0] << " " << mem.p_mess->potC[holy_grail][1] << endl;
		  }
	      }
	  }
	FileFractal << "calling power_spectrum from initial_forces " << length << endl;
	FileFractal << "sizes a " << variance_rho.size() << " " << variance_pot.size() << " " << variance_force.size() << " " << variance_force_s.size() << endl;
	power_spectrum(mem.p_mess->potC,length,variance_rho,variance_pot,variance_force,variance_force_s,lev,frac.get_density_0(),true,mem);
	//
	FileFractal << "back from power " << lev << endl;
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
		FileFractal << "dead parrot " << dead_parrot << endl;
	      }
	    else if(mem.norm_what == 1 || mem.norm_what==2)
	      {
		assert(0);
		var_obs_0=((1.0-a)*variance_force_s[n1]+a*variance_force_s[n1+1]);
	      }
	    else if(mem.norm_what == 3)
	      {
		assert(0);
		var_obs_0=((1.0-a)*variance_pot[n1]+a*variance_pot[n1+1]);
	      }
	    else
	      assert(0);
	  }
	//
	for(int i=0;i < length/2;i++)
	  {
	    FileFractal << "real variance a " << i << " " << (double)i/(double)length << " " << variance_rho[i] << endl;
	  }
	for(ptrdiff_t kx=mem.p_mess->start_x;kx < mem.p_mess->start_x+mem.p_mess->length_x ; kx++)
	  {
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
	  }
	mem.p_mess->fftw_complex_to_real();
	info_to_slices(mem,frac,lev);
	if(!lev==0)
	  {
	    for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
		group_itr!=mem.all_groups[lev].end();group_itr++)
	      {
		Group& group=**group_itr;
		potential_start(group);
	      }
	  }
	slices_to_pot_init(mem,frac,lev);
	//
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
	    group_itr!=mem.all_groups[lev].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    force_at_point(group,frac);
	  }
      }
    if(highest_level_fft < highest_level_used)
      {
	for(int lev=highest_level_fft+1;lev <= highest_level_used;++lev)
	  {
	    for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
		group_itr!=mem.all_groups[lev].end();group_itr++)
	      {
		Group& group=**group_itr;
		potential_start(group);
		force_at_point(group,frac);
	      } 
	  }
      }
    for(int level=0;level <= frac.get_level_max();level++)
      {
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	    group_itr!=mem.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    group.scale_pot_forces(dead_parrot);
	    double varx,vary,varz;
	    group.get_force_variance(varx,vary,varz);
	    FileFractal << "point variances " << group.get_level() << " " << group.list_points.size() << " " << varx << " " << vary << " " << varz << endl;
	    FileFractal << "sharpiea " << frac.get_level_max() << endl;
	    force_at_particle_sharp(group,frac); 		
	    FileFractal << "sharpieb " << frac.get_level_max() << endl;
	  }
      }
    FileFractal << "calling fractal " << &frac << endl;
    sum_pot_forces(frac);
    //  assert(0);
    FileFractal << "whatt0 " << endl;
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
	FileFractal << " level forces " << lev  << " " << sum1[0] << " " << sum2[0] << " " << sum1[1] << " " << sum2[1] << " " << sum1[2] << " " << sum2[2] << endl;
      } 
    sum_pot_forces(frac);
    //
    mem.p_mess->free_potRC();
    mem.p_mess->return_Slice_pos.clear();
    mem.p_mess->return_group.clear();
    mem.p_mess->return_point.clear();
    mem.p_mess->return_node.clear();
  }
}
