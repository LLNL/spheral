#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void shrink_cube(double SHRINK,vector <double>& xmin,vector <double>& xmax,Fractal_Memory* PFM,
		   vector <double>& posx,vector <double>& posy,vector <double>& posz,
		   int number_particles,vector <double>& xmini,vector <double>& xmaxy)
  {
    if(SHRINK < 1.0e-5)
      return;
    SHRINK=min(SHRINK,1.0);
    vector <double>minimax(6);
    auto mmx=std::minmax_element(posx.begin(),posx.begin()+number_particles);
    auto mmy=std::minmax_element(posy.begin(),posy.begin()+number_particles);
    auto mmz=std::minmax_element(posz.begin(),posz.begin()+number_particles);
    minimax[0]=-*mmx.first;
    minimax[1]=-*mmy.first;
    minimax[2]=-*mmz.first;
    minimax[3]=*mmx.second;
    minimax[4]=*mmy.second;
    minimax[5]=*mmz.second;
    PFM->p_mess->Find_Max_DOUBLE(minimax,6);
    minimax[0]*=-1.0;
    minimax[1]*=-1.0;
    minimax[2]*=-1.0;
    vector <double> dmaxy(3);
    vector <double> center(3);
    center[0]=(minimax[3]+minimax[0])*0.5;
    center[1]=(minimax[4]+minimax[1])*0.5;
    center[2]=(minimax[5]+minimax[2])*0.5;
    double dmax=max(max(minimax[3]-minimax[0],minimax[4]-minimax[1]),minimax[5]-minimax[2]);
    if(dmax >=xmax[0]-xmin[0])
      {
	xmini=xmin;
	xmaxy=xmax;
	return;
      }
    dmax=dmax*SHRINK+(xmax[0]-xmin[0])*(1.0-SHRINK);
    xmini[0]=center[0]-dmax*0.51;
    xmini[1]=center[1]-dmax*0.51;
    xmini[2]=center[2]-dmax*0.51;
    xmaxy[0]=center[0]+dmax*0.51;
    xmaxy[1]=center[1]+dmax*0.51;
    xmaxy[2]=center[2]+dmax*0.51;
    if(PFM->p_mess->FractalRank == 0)
      cerr << " Shrink Box " << PFM->steps << " " << xmini[0] << " " << xmini[1] << " " << xmini[2] << " "  << xmaxy[0] << " " << xmaxy[1] << " " << xmaxy[2] << " " << pow(dmax*1.02,3) << "\n";
  }
}
namespace FractalSpace
{
  void shrink_cube(Fractal_Memory* PFM,double SHRINK,vector <double>& xmin,vector <double>& xmax,
		   vector <double>& xmini,vector <double>& xmaxy)
  {
    if(SHRINK < 1.0e-5)
      return;
    SHRINK=min(SHRINK,1.0);
    Fractal* PF=PFM->p_fractal;
    vector <double>minimax(6,DBL_MIN);
    vector<double>pos(3);
    for(int ni=0;ni<PF->get_number_particles();ni++)
      {
	PF->particle_list[ni]->get_pos(pos);
	for(int k : {0,1,2})
	  {
	    minimax[k]=max(minimax[k],-pos[k]);
	    minimax[k+3]=max(minimax[k+3],pos[k]);
	  }
      }
    PFM->p_mess->Find_Max_DOUBLE(minimax,6);
    minimax[0]*=-1.0;
    minimax[1]*=-1.0;
    minimax[2]*=-1.0;
    vector <double> dmaxy(3);
    vector <double> center(3);
    center[0]=(minimax[3]+minimax[0])*0.5;
    center[1]=(minimax[4]+minimax[1])*0.5;
    center[2]=(minimax[5]+minimax[2])*0.5;
    double dmax=max(max(minimax[3]-minimax[0],minimax[4]-minimax[1]),minimax[5]-minimax[2]);
    if(dmax >=xmax[0]-xmin[0])
      {
	xmini=xmin;
	xmaxy=xmax;
	return;
      }
    dmax=dmax*SHRINK+(xmax[0]-xmin[0])*(1.0-SHRINK);
    xmini[0]=center[0]-dmax*0.51;
    xmini[1]=center[1]-dmax*0.51;
    xmini[2]=center[2]-dmax*0.51;
    xmaxy[0]=center[0]+dmax*0.51;
    xmaxy[1]=center[1]+dmax*0.51;
    xmaxy[2]=center[2]+dmax*0.51;
    if(PFM->p_mess->FractalRank == 0)
      cerr << " Shrink Box " << PFM->steps << " " << xmini[0] << " " << xmini[1] << " " << xmini[2] << " "  << xmaxy[0] << " " << xmaxy[1] << " " << xmaxy[2] << " " << pow(dmax*1.02,3) << "\n";
  }
}
