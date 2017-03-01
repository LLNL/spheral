#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void halo_force(Fractal& frac)
  {
    frac.copy_particle_list_to_list_sorted();
    sort(frac.particle_list_sorted.begin(),frac.particle_list_sorted.end(),rad_compare);
    int np=frac.get_number_particles();
    Particle* plast=frac.particle_list_sorted[np-1];
    //
    double moat=0.0;
    if(!frac.get_periodic()) moat=max(max(1,frac.get_moat_0()),frac.get_padding());
    double radmin=0.5-moat/(double)(frac.get_grid_length());
    if(plast->get_r2(0.5,0.5,0.5) <= radmin) return;
    //
    const double con10(sqrt(3.0));
    const double con11(sqrt(1.5));
    const double con20(sqrt(1.25));
    const double con21(sqrt(7.5));
    const double con22(sqrt(1.875));
    double A00(0.0);
    double A10(0.0);
    complex <double> A11(0.0);
    double A20(0.0);
    complex <double> A21(0.0);
    complex <double> A22(0.0);
    for(int n=0;n<np;n++)
      {
	vector <double> pos(3);
	vector <double> field(4);
	double mass;
	Particle* p=frac.particle_list_sorted[n];
	p->get_posm(pos,mass);
	double dx=pos[0]-0.5;
	double dy=pos[1]-0.5;
	double dz=pos[2]-0.5;
	double R2=dx*dx+dy*dy+1.0e-30;
	double R=sqrt(R2);
	double cphi=dx/R;
	double sphi=dy/R;
	complex <double> Cxphi1(cphi,sphi);
	complex <double> Cxphi2(Cxphi1*Cxphi1);
	complex <double> dCxphi1(-sphi,cphi);
	complex <double> dCxphi2(2.0*dCxphi1*dCxphi1);
	double r2=R2+dz*dz;
	double r=sqrt(r2);
	double r3=r2*r;
	double r4=r2*r2;
	double ctheta=dz/r;
	double stheta=R/r;
	double P00=1.0;
	double P10=ctheta*con10;
	double dP10=-stheta*con10;
	double P11=-stheta*con11;
	double dP11=-ctheta*con11;
	double P20=(3.0*ctheta*ctheta-1.0)*con20;
	double dP20=-6.0*stheta*ctheta*con20;
	double P21=-stheta*ctheta*con21;
	double dP21=(-ctheta*ctheta+stheta*stheta)*con21;
	double P22=stheta*stheta*con22;
	double dP22=2.0*ctheta*stheta*con22;
	if(r > radmin)
	  {
	    field[0]=-
	      (A00/r+A10*P10/r2+A20*P20/r3+
	       2.0*real(A11*P11*Cxphi1/r2+A21*P21*Cxphi1/r3+
			A22*P22*Cxphi2/r3));
	    double frad=-
	      (A00/r2+2.0*A10*P10/r3+3.0*A20*P20/r4+
	       2.0*real(2.0*A11*P11*Cxphi1/r3+3.0*A21*P21*Cxphi1/r4+
			3.0*A22*P22*Cxphi2/r4));
	    double ftheta=
	      (A10*dP10/r2+A20*dP20/r3+
	       2.0*real(A11*dP11*Cxphi1/r2+A21*dP21*Cxphi1/r3+
			A22*dP22*Cxphi2/r3))/r;
	    double fphi=
	      (2.0*real(A11*P11*dCxphi1/r2+A21*P21*dCxphi1/r3+
			A22*P22*dCxphi2/r3))*stheta/r;
	    // rvec=(dx,dy,dz)/r
	    //thetavec=(ctheta*cphi,ctheta*sphi,-stheta)
	    //phivec=(-sphi,cphi,0)
	    field[1]=frad*dx/r+ftheta*ctheta*cphi-fphi*sphi;
	    field[2]=frad*dy/r+ftheta*ctheta*sphi+fphi*cphi;
	    field[3]=frad*dz/r-ftheta*stheta;
	    p->add_field_pf(field);
	  }
	A00+=mass*P00;
	A10+=mass*r*P10;
	A11+=mass*r*P11*conj(Cxphi1);
	A20+=mass*r2*P20;
	A21+=mass*r2*P21*conj(Cxphi1);
	A22+=mass*r2*P22*conj(Cxphi2);
      }
    double B00(0.0);
    double B10(0.0);
    complex <double> B11(0.0);
    double B20(0.0);
    complex <double> B21(0.0);
    complex <double> B22(0.0);
    for(int n=np-1;n>=0;n--)
      {
	vector <double> pos(3);
	vector <double> field(4);
	double mass;
	Particle* p=frac.particle_list_sorted[n];
	p->get_posm(pos,mass);
	double dx=pos[0]-0.5;
	double dy=pos[1]-0.5;
	double dz=pos[2]-0.5;
	double R2=dx*dx+dy*dy+1.0e-30;
	double R=sqrt(R2);
	double cphi=dx/R;
	double sphi=dy/R;
	complex <double> Cxphi1(cphi,sphi);
	complex <double> Cxphi2(Cxphi1*Cxphi1);
	complex <double> dCxphi1(-sphi,cphi);
	complex <double> dCxphi2(2.0*dCxphi1*dCxphi1);
	double r2=R2+dz*dz;
	double r=sqrt(r2);
	double r3=r2*r;
	double ctheta=dz/r;
	double stheta=R/r;
	double P00=1.0;
	double P10=ctheta*con10;
	double dP10=-stheta*con10;
	double P11=-stheta*con11;
	double dP11=-ctheta*con11;
	double P20=(3.0*ctheta*ctheta-1.0)*con20;
	double dP20=-6.0*stheta*ctheta*con20;
	double P21=-stheta*ctheta*con21;
	double dP21=(-ctheta*ctheta+stheta*stheta)*con21;
	double P22=stheta*stheta*con22;
	double dP22=2.0*ctheta*stheta*con22;
	field[0]=-
	  (B00+B10*P10*r+B20*P20*r2+
	   2.0*real(B11*P11*Cxphi1*r+B21*P21*Cxphi1*r2+
		    B22*P22*Cxphi2*r2));
	double frad=
	  (B10*P10+2.0*B20*P20*r+
	   2.0*real(B11*P11*Cxphi1+2.0*B21*P21*Cxphi1*r+
		    2.0*B22*P22*Cxphi2*r));
	double ftheta=
	  (B10*dP10*r+B20*dP20*r2+
	   2.0*real(B11*dP11*Cxphi1*r+B21*dP21*Cxphi1*r2+
		    B22*dP22*Cxphi2*r2));
	double fphi=
	  (2.0*real(B11*P11*dCxphi1*r+B21*P21*dCxphi1*r2+
		    B22*P22*dCxphi2*r2));
	// rvec=(dx,dy,dz)/r
	//thetavec=(ctheta*cphi,ctheta*sphi,-stheta)
	//phivec=(-sphi,cphi,0)
	field[1]=frad*dx/r+ftheta*ctheta*cphi-fphi*sphi;
	field[2]=frad*dy/r+ftheta*ctheta*sphi+fphi*cphi;
	field[3]=frad*dz/r-ftheta*stheta;
	p->add_field_pf(field);
	if(r > radmin)
	  {
	    B00+=mass*P00/r;
	    B10+=mass*P10/r2;
	    B11+=mass*P11*conj(Cxphi1)/r2;
	    B20+=mass*P20/r3;
	    B21+=mass*P21*conj(Cxphi1)/r3;
	    B22+=mass*P22*conj(Cxphi2)/r3;
	  }
      }
    frac.particle_list_sorted.clear();
  }
  bool rad_compare(Particle* par1,Particle* par2)
  {
    return par1->get_r2(0.5,0.5,0.5) <= par2->get_r2(0.5,0.5,0.5);
  }
}
