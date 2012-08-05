#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void power_spectrum(fftw_complex* rhoC,const int& length,vector <double>& variance_rho,vector <double>& variance_pot,
		      vector <double>& variance_force,vector <double>& variance_force_s,const int& lev,const double& d0,const bool& do_var,
		      Fractal_Memory& mem)
  {
    cout << "sizes b " << variance_rho.size() << " " << variance_pot.size() << " " << variance_force.size() << " " << variance_force_s.size() << endl;
    static ofstream FilePow;
    if(!FilePow.is_open())
      FilePow.open("power.d");
    FilePow.precision(4);
    static ofstream FileVar;
    if(!FileVar.is_open())
      FileVar.open("variance.d");
    FileVar.precision(4);
    cout << "enter power_spectrum " << length << endl;
    vector <double> green_2(length+1);
    vector <double> force_1(length+1);
    vector <double> force_2(length+1);
    //  return;
    double var_0=pow((double)(length*length*length)*d0,-2);
    double pi=4.0*atan(1.0);
    double twopi=8.0*atan(1.0);
    double force_const=4.0*pi/(double)(length*length);
    double d_step_wave=pow(2.0,lev);
    int i_step_wave=Misc::pow(2,lev);
    double spam_6=2.0/pow((double)(length),6);
    cout << "spam_6= " << spam_6 << endl;
    for (int k=0;k <length;++k)
      {
	double aa=pi*(double)(k)/(double)(length);
	force_1[k]=(double)(length)*sin(2.0*aa);
	force_2[k]=2.0*(double)(length)*sin(aa);
	green_2[k]=pow(2.0*sin(aa),2);
	//      cout << "initial scalings " << k << " " << force_1[k] << " " << force_2[k] << endl;
      }
    const int nyq=length/2;
    int length_c=nyq+1;
    double length_inv=twopi/(double)length;
    assert(nyq > 0);
    int vec_length=length*i_step_wave;
    vector  <double> power(vec_length,0.0);
    vector  <double> sum_0(vec_length,1.0e-30);
    vector  <double> sum_1(vec_length,0.0);
    //
    int counts=0;
    for(int kx=0;kx < length;++kx)
      {
	const int ka=min(kx,length-kx);
	for(int ky=0;ky < length;++ky)
	  {
	    const int kb=min(ky,length-ky);
	    for(int kz=0;kz < nyq+1;++kz)
	      {
		const int kc=min(kz,length-kz);
		double k=sqrt((double)(ka*ka+kb*kb+kc*kc));
		int n=(int)k;
		int n_s=n*i_step_wave;
		int holy_grail=mem.fftw_where(kx,ky,kz,length,length_c);
		double square=pow(rhoC[holy_grail][0],2)+pow(rhoC[holy_grail][1],2);
		/*
		cout << "pow a " << " ";
		cout << kx << " ";
		cout << ky << " ";
		cout << kz << " ";
		cout << holy_grail << " ";
		cout << rhoC[holy_grail][0] << " ";
		cout << rhoC[holy_grail][1] << " ";
		cout << endl ;
		*/
		//		double square=norm(rho(kx,ky,kz));
		//	      cout << "pow haha " << ka << " "<< kb << " "<< kc << " "<< k << " "<< n << " " << square << endl;
		if(square > 0.0)
		  {
		    sum_0[n_s]+=1.0;
		    sum_1[n_s]+=k;
		    power[n_s]+=square;
		    double spam=k*length_inv;
		    double pot=square*pow(force_const/(green_2[ka]+green_2[kb]+green_2[kc]+1.0e-30),2);
		    double f1=(pow(force_1[ka],2)+pow(force_1[kb],2)+pow(force_1[kc],2))*pot;
		    double f2=(pow(force_2[ka],2)+pow(force_2[kb],2)+pow(force_2[kc],2))*pot;
		    if(do_var)
		      {
			for(int nv=0;nv < nyq;nv++)
			  {
			    counts++;
			    double filter=var_0*Misc::square_filter((double)(nv)*spam*d_step_wave);
			    variance_rho[nv]+=square*filter;
			    variance_pot[nv]+=pot*filter;
			    variance_force[nv]+=f1*filter;
			    variance_force_s[nv]+=f2*filter;
			  }
		      }
		  }
	      }
	  }
      }
    cout << "count power " << counts << endl;
    cout << "var zero " << variance_rho[0] << endl;
    for(int n=0;n <=nyq*i_step_wave;n+=i_step_wave)
      {
	int nv=n/i_step_wave;
	sum_1[n]/=sum_0[n];
	power[n]/=sum_0[n];
	//	cout << "pow and var " << lev << " " << n << " " << power[n] << " " << variance_rho[n] << endl;
	FilePow << n << "\t " << scientific << sum_1[n] << "\t " << power[n] << endl;
	if(do_var) FileVar << scientific << (double)(nv)/(double)length << "\t " << variance_rho[nv] << "\t " << variance_pot[nv] << "\t " 
			   << variance_force[nv] << "\t " << variance_force_s[nv] << endl;
      }
    cout << "leaving power " << lev << endl;
    //    assert(0);
  }
}
