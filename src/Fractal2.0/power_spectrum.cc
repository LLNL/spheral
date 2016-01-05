#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void power_spectrum(fftw_complex* rhoC,int length,vector <double>& variance_rho,vector <double>& variance_pot,
		      vector <double>& variance_force,vector <double>& variance_force_s,int lev,double d0,bool do_var,
		      Fractal_Memory& mem)
  {
    ofstream& FilePow=mem.p_fractal->p_file->DUMPS;
    ofstream& FileVar=mem.p_fractal->p_file->DUMPS;
//     FilePow << "sizes b " << variance_rho.size() << " " << variance_pot.size() << " " << variance_force.size() << " " << variance_force_s.size() << "\n";
//     FilePow << "enter power_spectrum " << length << "\n";
    vector <double> green_2(length+1);
    vector <double> force_1(length+1);
    vector <double> force_2(length+1);
    do_var=do_var && lev == 0;
    double var_0=pow((double)(length*length*length)*d0,-2);
    double pi=4.0*atan(1.0);
    double twopi=8.0*atan(1.0);
    double force_const=4.0*pi/(double)(length*length);
    double d_step_wave=pow(2.0,lev);
    int i_step_wave=Misc::pow(2,lev);
    //    double spam_6=2.0/pow((double)(length),6);
//     FilePow << "spam_6= " << spam_6 << "\n";
    for (int k=0;k <length;++k)
      {
	double aa=pi*(double)(k)/(double)(length);
	force_1[k]=(double)(length)*sin(2.0*aa);
	force_2[k]=2.0*(double)(length)*sin(aa);
	green_2[k]=pow(2.0*sin(aa),2);
	//      FilePow << "initial scalings " << k << " " << force_1[k] << " " << force_2[k] << "\n";
      }
    vector <double>sum_sq(length,0.0);
    vector <double>sum_p(length,0.0);
    vector <double>sum_f1(length,0.0);
    vector <double>sum_f2(length,0.0);
    const int nyq=length/2;
    int length_c=nyq+1;
    double length_inv=twopi/(double)length;
    assert(nyq > 0);
    int vec_length=length*i_step_wave;
    vector  <double> power(vec_length,0.0);
    vector  <double> sum_0(vec_length,1.0e-30);
    vector  <double> sum_1(vec_length,0.0);
    //
    int begin_x=mem.p_mess->start_x;
    int end_x=begin_x+mem.p_mess->length_x;
    int counts=0;
    for(int kx=begin_x;kx < end_x;++kx)
      {
	const int ka=min(kx,length-kx);
	for(int ky=0;ky < length;++ky)
	  {
	    const int kb=min(ky,length-ky);
	    for(int kz=0;kz < nyq+1;++kz)
	      {
		const int kc=min(kz,length-kz);
		double k=sqrt((double)(ka*ka+kb*kb+kc*kc));
		int nk=k;
		int n=static_cast<int>(k);
		int n_s=n*i_step_wave;
		int holy_grail=mem.fftw_where(kx,ky,kz,length,length_c);
		double square=pow(rhoC[holy_grail][0],2)+pow(rhoC[holy_grail][1],2);
		counts++;
		sum_sq[nk]+=square;
		//		FilePow << "power " << kx << " " << ky << " " << kz << " " << square << "\n";
		if(square > 0.0)
		  {
		    sum_0[n_s]+=1.0;
		    sum_1[n_s]+=k;
		    power[n_s]+=square;
		    //		    double spam=k*length_inv;
		    double pot=square*pow(force_const/(green_2[ka]+green_2[kb]+green_2[kc]+1.0e-30),2);
		    double f1=(pow(force_1[ka],2)+pow(force_1[kb],2)+pow(force_1[kc],2))*pot;
		    double f2=(pow(force_2[ka],2)+pow(force_2[kb],2)+pow(force_2[kc],2))*pot;
		    sum_p[nk]+=pot;
		    sum_f1[nk]+=f1;
		    sum_f2[nk]+=f2;
		  }
	      }
	  }
      }
    if(do_var)
      {
	for(int nv=0;nv < nyq;nv++)
	  {
	    double anv=nv;
	    for(int k=0;k<nyq;k++)
	      {
		double spam=static_cast<double>(k)*length_inv;
		double filter=var_0*Misc::square_filter(anv*spam*d_step_wave);
		variance_rho[nv]+=sum_sq[k]*filter;
		variance_pot[nv]+=sum_p[k]*filter;
		variance_force[nv]+=sum_f1[k]*filter;
		variance_force_s[nv]+=sum_f2[k]*filter;
	      }
	  }
      }
    vector <double>ssp(vec_length*3);
    for(int ni=0;ni<vec_length;ni++)
      {
	ssp[ni]=sum_0[ni];
	ssp[ni+vec_length]=sum_1[ni];
	ssp[ni+2*vec_length]=power[ni];
      }
    mem.p_mess->Find_Sum_DOUBLE(ssp,3*vec_length);
    //    mem.p_mess->Find_Sum_DOUBLE(sum_0,vec_length);
    //    mem.p_mess->Find_Sum_DOUBLE(sum_1,vec_length);
    //    mem.p_mess->Find_Sum_DOUBLE(power,vec_length);
    int how_long=nyq*i_step_wave+1;
    if(do_var)
      {
	mem.p_mess->Find_Sum_DOUBLE(variance_rho,how_long);
	mem.p_mess->Find_Sum_DOUBLE(variance_pot,how_long);
	mem.p_mess->Find_Sum_DOUBLE(variance_force,how_long);
	mem.p_mess->Find_Sum_DOUBLE(variance_force_s,how_long);
      }
    //    FilePow << "count power " << counts << "\n";
    //    FilePow << "var zero " << variance_rho[0] << "\n";
    for(int n=0;n <=nyq*i_step_wave;n+=i_step_wave)
      {
	int nv=n/i_step_wave;
	sum_1[n]/=sum_0[n];
	power[n]/=sum_0[n];
	if(Mess::IAMROOT)
	  FilePow << n << "\t " << scientific << sum_1[n] << "\t " << power[n] << "\n";
	if(Mess::IAMROOT)
	  if(do_var) 
	    FileVar << scientific << (double)(nv)/(double)length << "\t " << variance_rho[nv] << "\t " << variance_pot[nv] << "\t " 
		    << variance_force[nv] << "\t " << variance_force_s[nv] << "\n";
      }
    //    FilePow << "leaving power " << lev << "\n";
  }
}
