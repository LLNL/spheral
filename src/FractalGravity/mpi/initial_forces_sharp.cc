#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void initial_forces_sharp(Fractal_Memory& fractal_memory,Fractal& fractal,Misc& misc)
  {
    cout << "enter initial_forces " << endl;
    int highest_level_used=0;
    for(int level=0;level <= fractal.get_level_max();level++)
      {
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    highest_level_used=max(highest_level_used,group.get_level());
	  }
      }
    //  assert(0);
    int length=fractal.get_grid_length();
    int zoom=Misc::pow(2,fractal.get_level_max());
    vector <int>Box(6);
    fractal.getBox(Box);
    vector <int>Box3(3);
    Box3[0]=Box[0]*zoom;
    Box3[1]=Box[2]*zoom;
    Box3[2]=Box[4]*zoom;
    assert(length >0);
    int nyq=length/2;
    int length_c=nyq+1;
    int real_length=nyq*Misc::pow(2,highest_level_used)+1;
    vector <double> variance_rho(real_length,0.0);
    vector <double> variance_pot(real_length,0.0);
    vector <double> variance_force(real_length,0.0);
    vector <double> variance_force_s(real_length,0.0);
    double rand_max=(double)RAND_MAX;

    size_t sizeR=sizeof(double);
    double* potR;
    fftw_complex* potC;
    potR=(double*) fftw_malloc(sizeR*2*length_c*length*length);
    potC=(fftw_complex*) potR;
    fftw_plan plan_cr=fftw_plan_dft_c2r_3d(length,length,length,potC,potR,FFTW_ESTIMATE);
    double pi=4.0*atan(1.0);
    double twopi=2.0*pi;
    double fourpi=4.0*pi;
    assert(fractal_memory.highest_level_init >=0 && fractal_memory.highest_level_init <= fractal.get_level_max());
    int highest_level_fft=min(fractal_memory.highest_level_init,highest_level_used);
    double var_initial=fractal_memory.sigma_initial*fractal_memory.sigma_initial;
    double dead_parrot=0.0;
    cout << "fractal_memory.scaling= " << fractal_memory.scaling << " " << var_initial << endl;
    vector <double> green_2(length+1);
    for (int k=0;k <length;++k)
      {
	double aa=pi*(double)(k)/(double)(length);
	green_2[k]=pow(2.0*sin(aa),2);
      }
    //
    for(int lev=0;lev <= highest_level_fft;++lev)
      {
	cout << "make power " << lev << " " << highest_level_fft << endl;
	double boost_power=pow(8.0,lev);
	double boost_scale=pow(2.0,lev);
	double grav_const=fourpi/boost_scale/boost_scale/(double)(length*length);
	int division=Misc::pow(2,fractal.get_level_max()-lev);
	int wrapping=length*division;
	int cut_lev=-1;
	if(lev > 0) cut_lev=nyq/2;
	double step_wave=pow(2.0,lev);
	for(int kx=0;kx < length ; kx++)
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
		      amplitude_raw=sqrt(boost_power*cosmos_power(k/fractal_memory.scaling,fractal_memory));
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
		    double norwegian_blue=Fractal::my_rand_not_zero(rand_max);
		    double amplitude_random=amplitude_raw*sqrt(-2.0*log(norwegian_blue));
		    double pot_k=amplitude_random*step_wave;
		    int holy_grail=2*fftw_where(kx,ky,kz,length,length_c);
		    potR[holy_grail]=pot_k*cos(angle);
		    potR[holy_grail+1]=pot_k*cos(angle);
		    //		    pot(kx,ky,kz)=pot_k*Complex(cos(angle),sin(angle));
		  }
	      }
	  }
	cout << "calling power_spectrum from initial_forces " << length << endl;
	cout << "sizes a " << variance_rho.size() << " " << variance_pot.size() << " " << variance_force.size() << " " << variance_force_s.size() << endl;
	power_spectrum(potR,length,variance_rho,variance_pot,variance_force,variance_force_s,lev,fractal.get_density_0(),true);
	//
	cout << "back from power " << lev << endl;
	if(lev == 0)
	  {
	    double a=fractal_memory.norm_scale*(double)length;
	    int n1=(int)a;
	    a-=(double)n1;
	    double var_obs_0=-1.0;
	    if(fractal_memory.norm_what == 0)
	      {
		double var_obs=(1.0-a)*variance_rho[n1]+a*variance_rho[n1+1];
		dead_parrot=sqrt(var_initial/var_obs);
		cout << "dead parrot " << dead_parrot << endl;
	      }
	    else if(fractal_memory.norm_what == 1 || fractal_memory.norm_what==2)
	      {
		assert(0);
		var_obs_0=((1.0-a)*variance_force_s[n1]+a*variance_force_s[n1+1]);
	      }
	    else if(fractal_memory.norm_what == 3)
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
	    cout << "real variance a " << i << " " << (double)i/(double)length << " " << variance_rho[i] << endl;
	  }
	for(int kx=0;kx < length ; kx++)
	  {
	    int ka=min(kx,length-kx);
	    for(int ky=0;ky < length ; ky++)
	      {
		int kb=min(ky,length-ky);
		for(int kz=0;kz <= nyq ; kz++)
		  {
		    int kc=kz;
		    double g2=grav_const/(green_2[ka]+green_2[kb]+green_2[kc]+1.0e-30);
		    int brian=2*fftw_where(kx,ky,kz,length,length_c);
		    potR[brian]*=g2;
		    potR[brian+1]*=g2;
		    //		    double g2=green_2[ka]+green_2[kb]+green_2[kc]+1.0e-30;
		    //		    potC[fftw_where(kx,ky,kz,length,length_c)]*=grav_const/g2;
		  }
	      }
	  }
	//      pot(0,0,0)=1.0;
	fftw_execute(plan_cr);
	//	Backward3.fftNormalized(pot);
	//	cout << "real pot " << real(pot(0,0,0)) << " " << imag(pot(0,0,0))  << " " << real(pot(1,0,0)) << " " << imag(pot(1,0,0)) << endl;
	//      assert(0);
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[lev].begin();
	    group_itr!=fractal_memory.all_groups[lev].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    bool top_group= &group == misc.p_group_0;
	    cout << "group " << top_group << " " << &group << endl;
	    if(!top_group)
	      potential_start(group);
	    double sumpot=0.0;
	    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
	      {
		Point& point=**point_itr;
		int p_x=point.get_pos_point_x()-Box3[0];
		int p_y=point.get_pos_point_y()-Box3[1];
		int p_z=point.get_pos_point_z()-Box3[2];
		int p_xi=(p_x % wrapping)/division;
		int p_yi=(p_y % wrapping)/division;
		int p_zi=(p_z % wrapping)/division;
		if(top_group)
		  point.set_potential_point(potR[fftw_where(p_xi,p_yi,p_zi,length,length)]);
		else
		  point.add_potential_point(potR[fftw_where(p_xi,p_yi,p_zi,length,length)]);
		sumpot+=pow(point.get_potential_point(),2);
	      }
	    cout << "sumpot " << sumpot << endl;
	    force_at_point(group,fractal);
	  }
      }
    if(highest_level_fft < highest_level_used)
      {
	for(int lev=highest_level_fft+1;lev <= highest_level_used;++lev)
	  {
	    for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[lev].begin();
		group_itr!=fractal_memory.all_groups[lev].end();group_itr++)
	      {
		Group& group=**group_itr;
		if(lev == group.get_level())
		  {		  
		    potential_start(group);
		    force_at_point(group,fractal);
		  }
	      } 
	  }
      }
    for(int level=0;level <= fractal.get_level_max();level++)
      {
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    group.scale_pot_forces(dead_parrot);
	    double varx,vary,varz;
	    group.get_force_variance(varx,vary,varz);
	    cout << "point variances " << group.get_level() << " " << group.list_points.size() << " " << varx << " " << vary << " " << varz << endl;
	    cout << "sharpiea " << fractal.get_level_max() << endl;
	    force_at_particle_sharp(group,fractal); 		
	    cout << "sharpieb " << fractal.get_level_max() << endl;
	  }
      }
    cout << "calling fractal " << &fractal << endl;
    sum_pot_forces(fractal);
    //  assert(0);
    cout << "whatt0 " << endl;
    for(int lev=0;lev<=highest_level_used;lev++)
      {
	double sum0=1.0e-10;
	vector <double>sum1(3,1.0e-10);
	vector <double>sum2(3,1.0e-10);
	for(int i=0;i<fractal.get_number_particles();i++)
	  {
	    Particle& p=*fractal.particle_list[i];
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
	cout << " level forces " << lev  << " " << sum1[0] << " " << sum2[0] << " " << sum1[1] << " " << sum2[1] << " " << sum1[2] << " " << sum2[2] << endl;
      } 
    sum_pot_forces(fractal);
    //
    fftw_free(potR);
    potR=0;
    potC=0;
  }
  void sum_pot_forces(Fractal& fractal)
  {
    double n=fractal.get_number_particles();
    vector <double> sum_8(8,0.0);
    for(int p=0;p < n;++p)
      {
	Particle& particle=*fractal.particle_list[p];
	vector <double> field(4);
	particle.get_field_pf(field);
	sum_8[0]+=field[0];
	sum_8[1]+=pow(field[0],2);
	sum_8[2]+=field[1];
	sum_8[3]+=pow(field[1],2);
	sum_8[4]+=field[2];
	sum_8[5]+=pow(field[2],2);
	sum_8[6]+=field[3];
	sum_8[7]+=pow(field[3],2);
      }    
    for(int i=0;i < 4;i++)
      {
	sum_8[i*2]/=n;
	sum_8[i*2+1]/=n;
	sum_8[i*2+1]=sqrt(sum_8[i*2+1]-sum_8[i*2]*sum_8[i*2]);
      }
    cout << "sum_start " << sum_8[0] << " " << sum_8[1] << " " << sum_8[2] << " " << sum_8[3] << " " << sum_8[4] << " " << sum_8[5] << " " << sum_8[6] << " " << sum_8[7] << endl;
  }
}
