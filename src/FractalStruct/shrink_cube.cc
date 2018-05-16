#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  typedef deque<double>::iterator _ITD__;
  template <class ForwardIterator>
  void shrink_cube(double SHRINK,
		   const vector <double>& xmin,const vector <double>& xmax,
		   Fractal_Memory* PFM,
		   ForwardIterator posxb,
		   ForwardIterator posyb,
		   ForwardIterator poszb,
		   int number_particles,
		   vector <double>& xmini,vector <double>& xmaxy)
  {
    xmini=xmin;
    xmaxy=xmax;
    if(SHRINK < 1.0e-5)
      return;
    SHRINK=min(SHRINK,1.0);
    vector <double>minimax(6);
    vector<pair<int,int>>mm;
    auto result=std::minmax_element(posxb,posxb+number_particles);
    minimax[0]=-*result.first;
    minimax[3]=*result.second;
    result=std::minmax_element(posyb,posyb+number_particles);
    minimax[1]=-*result.first;
    minimax[4]=*result.second;
    result=std::minmax_element(poszb,poszb+number_particles);
    minimax[2]=-*result.first;
    minimax[5]=*result.second;
    PFM->p_mess->Find_Max_DOUBLE(minimax,6);
    double dmax=-1.0;
    vector <double> center(3);
    for(int ni:{0,1,2})
      {
	minimax[ni]*=-1.0;
	center[ni]=(minimax[ni+3]+minimax[ni])*0.5;
	dmax=max(dmax,minimax[ni+3]-minimax[ni]);
      }
    // if(dmax >=xmax[0]-xmin[0])
    //   {
    // 	xmini=xmin;
    // 	xmaxy=xmax;
    // 	return;
    //   }
    dmax=dmax*SHRINK+(xmax[0]-xmin[0])*(1.0-SHRINK);
    dmax=min(dmax,xmax[0]-xmin[0]);
    for(int ni : {0,1,2})
      {
	xmini[ni]=center[ni]-dmax*0.51;
	xmaxy[ni]=center[ni]+dmax*0.51;
      }
    if(PFM->p_mess->FractalRank == 0)
      cerr << " Shrink Cube " << PFM->steps << " " << xmini[0] << " " << xmini[1] << " " << xmini[2] << " "  << xmaxy[0] << " " << xmaxy[1] << " " << xmaxy[2] << " " << pow(dmax*1.02,3) << "\n";
  }
}
namespace FractalSpace
{
  template
  void shrink_cube(double SHRINK,
		   const vector <double>& xmin,const vector <double>& xmax,
		   Fractal_Memory* PFM,
		   _ITD__ posxb,
		   _ITD__ posyb,
		   _ITD__ poszb,
		   int number_particles,
		   vector <double>& xmini,vector <double>& xmaxy);
}
namespace FractalSpace
{
  void shrink_cube(double SHRINK,
		   const vector <double>& xmin,const vector <double>& xmax,Fractal_Memory* PFM,
		   vector <double>& posx,vector <double>& posy,vector <double>& posz,
		   int number_particles,vector <double>& xmini,vector <double>& xmaxy)
  {
    xmini=xmin;
    xmaxy=xmax;
    if(SHRINK < 1.0e-5)
      return;
    SHRINK=min(SHRINK,1.0);
    vector <double>minimax(6);
    vector<pair<int,int>>mm;
    auto result=std::minmax_element(posx.begin(),posx.begin()+number_particles);
    minimax[0]=-*result.first;
    minimax[3]=*result.second;
    result=std::minmax_element(posy.begin(),posy.begin()+number_particles);
    minimax[1]=-*result.first;
    minimax[4]=*result.second;
    result=std::minmax_element(posz.begin(),posz.begin()+number_particles);
    minimax[2]=-*result.first;
    minimax[5]=*result.second;
    PFM->p_mess->Find_Max_DOUBLE(minimax,6);
    double dmax=-1.0;
    vector <double> center(3);
    for(int ni:{0,1,2})
      {
	minimax[ni]*=-1.0;
	center[ni]=(minimax[ni+3]+minimax[ni])*0.5;
	dmax=max(dmax,minimax[ni+3]-minimax[ni]);
      }
    if(dmax >=xmax[0]-xmin[0])
      {
	xmini=xmin;
	xmaxy=xmax;
	return;
      }
    dmax=dmax*SHRINK+(xmax[0]-xmin[0])*(1.0-SHRINK);
    for(int ni : {0,1,2})
      {
	xmini[ni]=center[ni]-dmax*0.51;
	xmaxy[ni]=center[ni]+dmax*0.51;
      }
    if(PFM->p_mess->FractalRank == 0)
      cerr << " Shrink Box " << PFM->steps << " " << xmini[0] << " " << xmini[1] << " " << xmini[2] << " "  << xmaxy[0] << " " << xmaxy[1] << " " << xmaxy[2] << " " << pow(dmax*1.02,3) << "\n";
  }
}
namespace FractalSpace
{
  void FractalCube(Fractal_Memory* PFM,double SHRINK,
		   const vector <double>& xmin,const vector <double>& xmax,
		   vector <double>& xmini,vector <double>& xmaxy)
  {
    xmini=xmin;
    xmaxy=xmax;
    double dmax=pow(xmaxy[0]-xmini[0],3);
    Fractal* PF=PFM->p_fractal;
    if(SHRINK > 1.0e-5)
      {
	SHRINK=min(SHRINK,1.0);
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
	vector <double> center(3);
	dmax=-1.0;
	for(int ni : {0,1,2})
	  {
	    minimax[ni]*=-1.0;
	    center[ni]=(minimax[ni]+minimax[ni+3])*0.5;
	    dmax=max(dmax,minimax[ni+3]-minimax[ni]);
	  }
	if(dmax >=xmax[0]-xmin[0])
	  {
	    xmini=xmin;
	    xmaxy=xmax;
	  }
	else
	  {
	    dmax=dmax*SHRINK+(xmax[0]-xmin[0])*(1.0-SHRINK);
	    for(int ni : {0,1,2})
	      {
		xmini[ni]=center[ni]-dmax*0.51;
		xmaxy[ni]=center[ni]+dmax*0.51;
	      }
	  }
      }
    double dinv=1.0/(xmaxy[0]-xmini[0]);
    for(int ni=0;ni<PF->get_number_particles();ni++)
      {
	vector<double>pos(3);
	PF->particle_list[ni]->get_pos(pos);
	for(int k : {0,1,2})
	  pos[k]=(pos[k]-xmini[k])*dinv;
	PF->particle_list[ni]->set_pos(pos);
      }
    if(PFM->p_mess->FractalRank == 0)
      cerr << " Fractal Cube " << PFM->steps << " " << xmini[0] << " " << xmini[1] << " " << xmini[2] << " "  << xmaxy[0] << " " << xmaxy[1] << " " << xmaxy[2] << " " << pow(dmax*1.02,3) << "\n";
  }
}
