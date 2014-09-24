#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void shrink_cube(double SHRINK,vector <double>& xmin,vector <double>& xmax,Fractal_Memory* PFM,
		   vector <double>& posx,vector <double>& posy,vector <double>& posz,
		   int number_particles,vector <double>& xmini,vector <double>& xmaxy)
  {
    SHRINK=max(SHRINK,0.0);
    SHRINK=min(SHRINK,1.0);
    if(SHRINK < 1.0e-5)
      return;
    vector <double>minimax(6,-DBL_MAX);
    for(int ni=0;ni<number_particles;ni++)
      {
	minimax[0]=max(minimax[0],-posx[ni]);
	minimax[1]=max(minimax[1],-posy[ni]);
	minimax[2]=max(minimax[2],-posz[ni]);
	minimax[3]=max(minimax[3],posx[ni]);
	minimax[4]=max(minimax[4],posy[ni]);
	minimax[5]=max(minimax[5],posz[ni]);
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
    double dmax=minimax[3]-minimax[0];
    dmax=max(minimax[4]-minimax[1],dmax);
    dmax=max(minimax[5]-minimax[2],dmax);
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
      cout << " Shrink Box " << xmini[0] << " " << xmini[1] << " " << xmini[2] << " "  << xmaxy[0] << " " << xmaxy[1] << " " << xmaxy[2] << " " << pow(dmax*1.02,3) << "\n";
  }
}
