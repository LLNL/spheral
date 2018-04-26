#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  typedef deque<double>::iterator _ITD__;
  template <class ForwardIterator>
  void am_I_conservative_enough_isol(Fractal_Memory* PFM,ForwardIterator massesb,double G,
				     vector <double>& xmin,vector <double>& xmax,double correction,
				     ForwardIterator posxb,ForwardIterator posyb,
				     ForwardIterator poszb,ForwardIterator velxb,
				     ForwardIterator velyb,ForwardIterator velzb)
  {
    ofstream& FileEnergy=PFM->p_file->FileEnergy;
    ofstream& FileMom=PFM->p_file->DUMPS;
    int stride=100;
    int NP=PFM->number_particles;
    double dt5=correction*PFM->step_length;
    double vx,vy,vz;
    double pe=0.0;
    double ke=0.0;
    double p0=0.0;
    double p1=0.0;
    double p2=0.0;
    double m0=0.0;
    double m1=0.0;
    double m2=0.0;
    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	vector <double>pot(many);
	vector <double>fx(many);
	vector <double>fy(many);
	vector <double>fz(many);
	get_field(PFM,ni,many,G,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    if(!I_am_a_real_particle(PFM,nip))
	       continue;
	    vx=*velxb+fx[p]*dt5;
	    vy=*velyb+fy[p]*dt5;
	    vz=*velzb+fz[p]*dt5;
	    pe+=*massesb*pot[p];
	    ke+=*massesb*(vx*vx+vy*vy+vz*vz);
	    p0+=*massesb*vx;
	    p1+=*massesb*vy;
	    p2+=*massesb*vz;
	    m0+=*massesb*(*posyb*vz-*poszb*vy);
	    m1+=*massesb*(*poszb*vx-*posxb*vz);
	    m2+=*massesb*(*posxb*vy-*posyb*vx);
	    posxb++;
	    posyb++;
	    poszb++;
	    velxb++;
	    velyb++;
	    velzb++;
	    massesb++;
	  }
      }
    pe*=0.5;
    ke*=0.5;
    double te=pe+ke;
    vector <double>sums{te,pe,ke,p0,p1,p2,m0,m1,m2};
    PFM->p_mess->Find_Sum_DOUBLE(sums,9);
    FileEnergy << scientific << PFM->time <<  "\t " << PFM->steps << "\t " << 
      te << "\t " << pe << "\t " << ke << "\t " << sums[0] << "\t" << sums[1] << "\t" << sums[2] << "\n";
    FileMom << scientific << PFM->time << "\t " << PFM->steps << "\t " << 
      p0 << "\t " << p1 << "\t " << p2 << "\t " <<
      m0 << "\t " << m1 << "\t " << m2 << "\t " << 
      sums[3] << "\t" << sums[4] << "\t" << sums[5] << "\t" <<
      sums[6] << "\t" << sums[7] << "\t" << sums[8] << "\n";
  }
}
namespace FractalSpace
{
  template
  void am_I_conservative_enough_isol(Fractal_Memory* PFM,_ITD__ massesb,double G,
				     vector <double>& xmin,vector <double>& xmax,double correction,
				     _ITD__ posxb,_ITD__ posyb,
				     _ITD__ poszb,_ITD__ velxb,
				     _ITD__ velyb,_ITD__ velzb);
}
namespace FractalSpace
{
  void am_I_conservative_enough_isol(Fractal_Memory* PFM,vector <double>& masses,double G,
				     vector <double>& xmin,vector <double>& xmax,double correction,
				     vector <double>& posx,vector <double>& posy,vector <double>& posz,
				     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    ofstream& FileEnergy=PFM->p_file->FileEnergy;
    ofstream& FileMom=PFM->p_file->DUMPS;
    int stride=100;
    int NP=PFM->number_particles;
    double dt5=correction*PFM->step_length;
    vector <double>pot(stride);
    vector <double>fx(stride);
    vector <double>fy(stride);
    vector <double>fz(stride);
    double vx,vy,vz;
    double pe=0.0;
    double ke=0.0;
    double p0=0.0;
    double p1=0.0;
    double p2=0.0;
    double m0=0.0;
    double m1=0.0;
    double m2=0.0;
    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	get_field(PFM,ni,many,G,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    if(!I_am_a_real_particle(PFM,nip))
	       continue;
	    vx=velx[nip]+fx[p]*dt5;
	    vy=vely[nip]+fy[p]*dt5;
	    vz=velz[nip]+fz[p]*dt5;
	    pe+=masses[nip]*pot[p];
	    ke+=masses[nip]*(vx*vx+vy*vy+vz*vz);
	    p0+=masses[nip]*vx;
	    p1+=masses[nip]*vy;
	    p2+=masses[nip]*vz;
	    m0+=masses[nip]*(posy[nip]*vz-posz[nip]*vy);
	    m1+=masses[nip]*(posz[nip]*vx-posx[nip]*vz);
	    m2+=masses[nip]*(posx[nip]*vy-posy[nip]*vx);
	  }
      }
    pe*=0.5;
    ke*=0.5;
    double te=pe+ke;
    vector <double>sums{te,pe,ke,p0,p1,p2,m0,m1,m2};
    PFM->p_mess->Find_Sum_DOUBLE(sums,9);
    FileEnergy << scientific << PFM->time <<  "\t " << PFM->steps << "\t " << 
      te << "\t " << pe << "\t " << ke << "\t " << sums[0] << "\t" << sums[1] << "\t" << sums[2] << "\n";
    FileMom << scientific << PFM->time << "\t " << PFM->steps << "\t " << 
      p0 << "\t " << p1 << "\t " << p2 << "\t " <<
      m0 << "\t " << m1 << "\t " << m2 << "\t " << 
      sums[3] << "\t" << sums[4] << "\t" << sums[5] << "\t" <<
      sums[6] << "\t" << sums[7] << "\t" << sums[8] << "\n";
  }
}
