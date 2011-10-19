 #include "libs.hh"
#include "classes.hh"
#include "headers.hh"
void initial_forces(Fractal& fractal,Misc& misc,list <Chain*>& list_chains)
{
  cout << "enter initial_forces " << endl;
  int highest_level_used=0;
  for(list <Chain*>::const_iterator chain_itr=list_chains.begin();chain_itr != list_chains.end();++chain_itr)
    {
      Chain& chain=**chain_itr;
      for(list <Group*>::const_iterator group_itr=chain.list_groups.begin();group_itr != chain.list_groups.end();++group_itr)
	{
	  Group& group=**group_itr;
	  highest_level_used=max(highest_level_used,group.get_level());
	}
    }
  //  assert(0);
  double rand_max=(double)RAND_MAX;
  double rand_min=1.0/rand_max;
  int length=fractal.get_grid_length();
  assert(length >0);
  int nyq=length/2;
  size_t align=sizeof(Complex);
  array3 <Complex> pot(length,length,nyq+1,align);
  crfft3d Backward3 (pot);
  double pi=4.0*atan(1.0);
  double twopi=2.0*pi;
  double fourpi=4.0*pi;
  vector <double> force_1(length+1);
  vector <double> green_2(length+1);
  vector <double> scaling_2(length+1);
  vector <double> scaling_4(length+1);
  for (int k=0;k <length;++k)
    {
      double aa=pi*(double)(k)/(double)(length);
      force_1[k]=((double)length)*sin(2.0*aa);
      green_2[k]=pow(2.0*sin(aa),2);
      scaling_2[k]=Misc::sinc_2(2.0*aa);
      scaling_4[k]=Misc::sinc_2(4.0*aa);
      cout << "initial scalings " << k << " " << scaling_2[k] << " " << scaling_4[k] << endl;
    }
  int highest_level_init=fractal.get_parameters_int(513);
  int norm_what=fractal.get_parameters_int(514);
  assert(highest_level_init >=0 && highest_level_init <= fractal.get_level_max());
  int highest_level_fft=min(highest_level_init,highest_level_used);
  double n3=(double)(length*length*length);
  double norm_scale=fractal.get_parameters_double(514);
  double var_initial=fractal.get_parameters_double(515);
  double gauss_filter=fractal.get_parameters_bool(514);
  double  sum_rho=0.0;
  double  var_rho=0.0;
  double  sum_pot=0.0;
  double var_pot=0.0;
  double sum_force_x=0.0;
  double var_force_x=0.0;
  double sum_force_y=0.0;
  double var_force_y=0.0;
  double sum_force_z=0.0;
  double var_force_z=0.0;
  double sum_rho_0=0.0;
  double var_rho_0=0.0;
  double sum_pot_0=0.0;
  double var_pot_0=0.0;
  double sum_force_x_0=0.0;
  double var_force_x_0=0.0;
  double sum_force_y_0=0.0;
  double var_force_y_0=0.0;
  double sum_force_z_0=0.0;
  double var_force_z_0=0.0;
  //
  for(int lev=0;lev <= highest_level_fft;++lev)
    {
      bool on_top= lev == 0;
      double boost_power=pow(8.0,lev);
      double boost_scale=pow(2.0,lev);
      bool grav_const=(double)(length*length)*fourpi*boost_scale;
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
		  double k=step_wave*(double)(sqrt(kx*kx+ky*ky+kz*kz));
		  bool shorty=(shorty_x && shorty_y && shorty_z) || k < 0.001;
		  double scaling_shorty=1.0;
		  if(shorty)
		    scaling_shorty=1.0-(scaling_4[ka]+scaling_4[kb]+scaling_4[kc])/(scaling_2[ka]+scaling_2[kb]+scaling_2[kc]);
		  double amplitude_raw=0.0;
		  if(k > 0.001)
		    amplitude_raw=sqrt(scaling_shorty*boost_power*cosmos_power(k,fractal));
		  double angle=0.0;
		  if(!nyquist)
		    angle=twopi*rand()/rand_max;
		  double norwegian_blue=max(rand_min,rand()/rand_max);
		  double amplitude_random=amplitude_raw*sqrt(-2.0*log(norwegian_blue));
		  double g2=green_2[ka]+green_2[kb]+green_2[kc]+1.0e-30;
		  double pot_k=grav_const*amplitude_random*step_wave/g2;
		  pot(kx,ky,kz)=Complex(cos(angle),sin(angle))*pot_k;
		  double rho_k=amplitude_random;
		  double force_x=force_1[ka]*pot_k*step_wave;
		  double force_y=force_1[kb]*pot_k*step_wave;
		  double force_z=force_1[kc]*pot_k*step_wave;
		  double q=k*norm_scale;
		  double filter=1.0;
		  if(gauss_filter)
		    filter=exp(-q*q*0.5);
		  else
		    if(q > 0.001) filter=3.0*(sin(q)-q*cos(q))/(q*q*q);
		  double filter_2=filter*filter;
		  //
		  //		  cout << kx << " "<< ky << " "<< kz << " " << rho_k << " " << pot_k << " " << scaling_shorty  <<  endl;
		  sum_rho+=rho_k*filter;
		  var_rho+=rho_k*rho_k*filter_2;
		  sum_pot+=pot_k*filter;
		  var_pot+=pot_k*pot_k*filter_2;
		  sum_force_x+=force_x*filter;
		  var_force_x+=force_x*force_x*filter_2;
		  sum_force_y+=force_y*filter;
		  var_force_y+=force_y*force_y*filter_2;
		  sum_force_z+=force_z*filter;
		  var_force_z+=force_z*force_z*filter*filter_2;
		  if(on_top)
		    {
		      sum_rho_0+=rho_k*filter;
		      var_rho_0+=rho_k*rho_k*filter_2;
		      sum_pot_0+=pot_k*filter;
		      var_pot_0+=pot_k*pot_k*filter_2;
		      sum_force_x_0+=force_x*filter;
		      var_force_x_0+=force_x*force_x*filter_2;
		      sum_force_y_0+=force_y*filter;
		      var_force_y_0+=force_y*force_y*filter_2;
		      sum_force_z_0+=force_z*filter;
		      var_force_z_0+=force_z*force_z*filter_2;
		    }
		}
	    }
	}
      cout << "calling power_spectrum from initial_forces " << length << endl;
      power_spectrum(pot,length);
      Backward3.fftNormalized(pot);
      //
      double conversion_f=(double)(length*Misc::pow(2,lev)/2);
      for(list <Chain*>::const_iterator chain_itr=list_chains.begin();chain_itr != list_chains.end();++chain_itr)
	{
	  Chain& chain=**chain_itr;
	  for(list <Group*>::const_iterator group_itr=chain.list_groups.begin();group_itr != chain.list_groups.end();++group_itr)
	    {
	      Group& group=**group_itr;
	      if(lev ==group.get_level())
		{
		  bool top_group= &group == misc.p_group_0;
		  if(!top_group)
		    {
		      potential_start(group,misc); 		
		      force_at_point_initial(group,fractal);
		    }
		  for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
		    {
		      Point& point=**point_itr;
		      int p_x=point.get_pos_point_x();
		      int p_y=point.get_pos_point_y();
		      int p_z=point.get_pos_point_z();
		      int p_xi=(p_x % wrapping)/division;
		      int p_xi_up=(p_xi+1) % length;
		      int p_xi_down=(p_xi-1+length) % length;
		      int p_yi=(p_y % wrapping)/division;
		      int p_yi_up=(p_yi+1) % length;
		      int p_yi_down=(p_yi-1+length) % length;
		      int p_zi=(p_z % wrapping)/division;
		      int first=p_zi/2;
		      if(top_group)
			{
			  if (p_zi % 2 == 0)
			    {
			      point.set_potential_point(real(pot(p_xi,p_yi,first)));
			      point.set_force_point_x(real(pot(p_xi_down,p_yi,first)-pot(p_xi_up,p_yi,first))*conversion_f);
			      point.set_force_point_y(real(pot(p_xi,p_yi_down,first)-pot(p_xi,p_yi_up,first))*conversion_f);
			      int first_down=(first-1+nyq) % nyq;
			      point.set_force_point_z(imag(pot(p_xi,p_yi,first_down)-pot(p_xi,p_yi,first))*conversion_f);
			    }
			  else
			    {
			      point.set_potential_point(imag(pot(p_xi,p_yi,first)));
			      point.set_force_point_x(imag(pot(p_xi_down,p_yi,first)-pot(p_xi_up,p_yi,first))*conversion_f);
			      point.set_force_point_y(imag(pot(p_xi,p_yi_down,first)-pot(p_xi,p_yi_up,first))*conversion_f);
			      int first_up=(first+1) % nyq;
			      point.set_force_point_z(real(pot(p_xi,p_yi,first)-pot(p_xi,p_yi,first_up))*conversion_f);
			    }
			}
		      else
			{
			  if (p_zi % 2 == 0)
			    {
			      point.add_potential_point(real(pot(p_xi,p_yi,first)));
			      point.add_force_point_x(real(pot(p_xi_down,p_yi,first)-pot(p_xi_up,p_yi,first))*conversion_f);
			      point.add_force_point_y(real(pot(p_xi,p_yi_down,first)-pot(p_xi,p_yi_up,first))*conversion_f);
			      int first_down=(first-1+nyq) % nyq;
			      point.add_force_point_z(imag(pot(p_xi,p_yi,first_down)-pot(p_xi,p_yi,first))*conversion_f);
			    }
			  else
			    {
			      point.add_potential_point(imag(pot(p_xi,p_yi,first)));
			      point.add_force_point_x(imag(pot(p_xi_down,p_yi,first)-pot(p_xi_up,p_yi,first))*conversion_f);
			      point.add_force_point_y(imag(pot(p_xi,p_yi_down,first)-pot(p_xi,p_yi_up,first))*conversion_f);
			      int first_up=(first+1) % nyq;
			      point.add_force_point_z(real(pot(p_xi,p_yi,first)-pot(p_xi,p_yi,first_up))*conversion_f);
			    }
			}
		    }
		}
	    }
	}      
    }
  if(highest_level_fft < highest_level_used)
    {
      for(int lev=highest_level_fft+1;lev <= highest_level_used;++lev)
	{
	  for(list <Chain*>::const_iterator chain_itr=list_chains.begin();chain_itr != list_chains.end();++chain_itr)
	    {
	      Chain& chain=**chain_itr;
	      for(list <Group*>::const_iterator group_itr=chain.list_groups.begin();group_itr != chain.list_groups.end();++group_itr)
		{
		  Group& group=**group_itr;
		  if(lev == group.get_level())
		    {
		      potential_start(group,misc); 		
		      force_at_point_initial(group,fractal);
		    }
		} 
	    }
	}
    }
  for(list <Chain*>::const_iterator chain_itr=list_chains.begin();chain_itr != list_chains.end();++chain_itr)
    {
      Chain& chain=**chain_itr;
      for(list <Group*>::const_iterator group_itr=chain.list_groups.begin();group_itr != chain.list_groups.end();++group_itr)
	{
	  Group& group=**group_itr;
	  force_at_particle(group,fractal,misc); 		
	}
    }
  fractal.sum_pot_forces();
  double var_obs_0=-1.0;
  if(norm_what==0)
    {
      sum_rho_0/=n3;
      var_obs_0=var_rho_0/n3-sum_rho_0*sum_rho_0;
    }
  else if(norm_what==1)
    {
      sum_force_x_0/=n3;
      sum_force_y_0/=n3;
      sum_force_z_0/=n3;
      var_obs_0=(var_force_x_0/n3-sum_force_x_0*sum_force_x_0+
	var_force_y_0/n3-sum_force_y_0*sum_force_y_0+
		 var_force_z_0/n3-sum_force_z_0*sum_force_z_0)/3.0;
    }
  else if(norm_what==2)
    {
      sum_pot_0/=n3;
      var_obs_0=var_pot_0/n3-sum_pot_0*sum_pot_0;
    }
  else
    {
      assert(0);
    }
  double dead_parrot=sqrt(var_initial/var_obs_0);;
  for(int p=0;p < fractal.get_number_particles();++p)
    {
      fractal.scale_pot_forces(p,dead_parrot);
    }
  fractal.sum_pot_forces();
}
