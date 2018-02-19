#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  MPI_Comm Fractal_Memory::FRACTAL_UNIVERSE;
  void Fractal_Memory::set_G(double Cavendish)
  {
    G=Cavendish;
  }
  void Fractal_Memory::get_G(double& Cavendish) const
  {
    Cavendish=G;
  }
  void Fractal_Memory::set_xmin(vector <double>& xm)
  {
    xmin=xm;
  }
  void Fractal_Memory::set_xmax(vector <double>& xm)
  {
    xmax=xm;
  }
  void Fractal_Memory::get_xmin(vector <double>& xm) const
  {
    xm=xmin;
  }
  void Fractal_Memory::get_xmax(vector <double>& xm) const
  {
    xm=xmax;
  }
  void Fractal_Memory::calc_FractalNodes()
  {
    FractalNodes=FractalNodes0*FractalNodes1*FractalNodes2;
    Boxes.resize(FractalNodes);
    int length=grid_length;
    if(!periodic)
      length++;
    int count=0;
    int j2b=-1;
    for(int m2=0;m2<FractalNodes2;m2++)
      {
	int j2a=j2b+1;
	j2b=((m2+1)*length)/FractalNodes2-1;
	int j1b=-1;
	for(int m1=0;m1<FractalNodes1;m1++)
	  {
	    int j1a=j1b+1;
	    j1b=((m1+1)*length)/FractalNodes1-1;
	    int j0b=-1;
	    for(int m0=0;m0<FractalNodes0;m0++)
	      {
		Boxes[count].resize(6);
		int j0a=j0b+1;
		j0b=((m0+1)*length)/FractalNodes0-1;
		Boxes[count][0]=j0a;
		Boxes[count][1]=j0b;
		Boxes[count][2]=j1a;
		Boxes[count][3]=j1b;
		Boxes[count][4]=j2a;
		Boxes[count][5]=j2b;
		count++;
	      }
	  }
      }
    int RANKY;
    MPI_Comm_rank(FractalWorld,&RANKY);
    for(int FR=0;FR < FractalNodes;FR++)
      {
	if(RANKY == 0)
	  p_file->FileFractal << "Box fracAA " << Boxes[FR][0] << " " << Boxes[FR][1] << " " << Boxes[FR][2] << " " 
			      << Boxes[FR][3] << " " << Boxes[FR][4] << " " << Boxes[FR][5] << "\n";

      }
    assert(count==FractalNodes);
  }
  void Fractal_Memory::calc_Buffers_and_more()
  {
//     vector < vector < vector <int> > > BoxesLev;
//     vector < vector < vector <int> > > BBoxesLev;
    Buffers.resize(FractalNodes);
    BBoxes.resize(FractalNodes);
    PBoxes.resize(FractalNodes);
    PBoxesLength.resize(FractalNodes);
    BoxesLev.resize(FractalNodes);
    BBoxesLev.resize(FractalNodes);
    PBoxesLev.resize(FractalNodes);
//     HRBoxesLev.resize(FractalNodes);
//     HSBoxesLev.resize(FractalNodes);

    for(int FR=0;FR<FractalNodes;FR++)
      {
	Buffers[FR].resize(6);
	BBoxes[FR].resize(6);
	PBoxes[FR].resize(6);
	PBoxesLength[FR].resize(3);
	for(int n=0;n<3;n++)
	  {
	    if(Boxes[FR][2*n] == 0 && !periodic)
	      Buffers[FR][2*n]=0;
	    else
	      Buffers[FR][2*n]=1;
	    if(Boxes[FR][2*n+1] == grid_length-1 && !periodic)
	      Buffers[FR][2*n+1]=0;
	    else
	      Buffers[FR][2*n+1]=1;
	    BBoxes[FR][2*n]=Boxes[FR][2*n]-Buffers[FR][2*n];
	    BBoxes[FR][2*n+1]=Boxes[FR][2*n+1]+Buffers[FR][2*n+1];
	    PBoxes[FR][2*n]=BBoxes[FR][2*n]-Buffers[FR][2*n];
	    PBoxes[FR][2*n+1]=BBoxes[FR][2*n+1]+Buffers[FR][2*n+1];
	    PBoxesLength[FR][n]=PBoxes[FR][2*n+1]-PBoxes[FR][2*n]+1;
	  }
	BoxesLev[FR].resize(level_max+1);
	BBoxesLev[FR].resize(level_max+1);
	PBoxesLev[FR].resize(level_max+1);
// 	HRBoxesLev[FR].resize(level_max+1);
// 	HSBoxesLev[FR].resize(level_max+1);
	BoxesLev[FR][0].resize(6);
	BBoxesLev[FR][0].resize(6);
	PBoxesLev[FR][0].resize(6);
	int zoom=Misc::pow(2,level_max);
	for(int n=0;n<6;n++)
	  {
	    BoxesLev[FR][0][n]=Boxes[FR][n]*zoom;
	    BBoxesLev[FR][0][n]=BBoxes[FR][n]*zoom;
	    PBoxesLev[FR][0][n]=PBoxes[FR][n]*zoom;
	  }
	for(int lev=1;lev<=level_max;lev++)
	  {
	    BoxesLev[FR][lev].resize(6);
	    BBoxesLev[FR][lev].resize(6);
	    PBoxesLev[FR][lev].resize(6);
	    PBoxesLev[FR][lev].resize(6);
// 	    HRBoxesLev[FR][lev].resize(6);
// 	    HSBoxesLev[FR][lev].resize(6);
	    zoom=Misc::pow(2,level_max-lev);
	    for(int n=0;n<3;n++)
	      {
		BoxesLev[FR][lev][2*n]=BoxesLev[FR][lev-1][2*n];
		BBoxesLev[FR][lev][2*n+1]=BBoxesLev[FR][lev-1][2*n+1];
		  
		BoxesLev[FR][lev][2*n+1]=BBoxesLev[FR][lev][2*n+1]-zoom*Buffers[FR][2*n+1];
		BBoxesLev[FR][lev][2*n]=BoxesLev[FR][lev][2*n]-zoom*Buffers[FR][2*n];
		  
		PBoxesLev[FR][lev][2*n+1]=BBoxesLev[FR][lev][2*n+1]+zoom*Buffers[FR][2*n+1];
		PBoxesLev[FR][lev][2*n]=BBoxesLev[FR][lev][2*n]-zoom*Buffers[FR][2*n];
		  
// 		HRBoxesLev[FR][lev][2*n+1]=PBoxesLev[FR][lev][2*n+1];
// 		HRBoxesLev[FR][lev][2*n]=BBoxesLev[FR][lev][2*n];
		  
// 		HSBoxesLev[FR][lev][2*n+1]=BoxesLev[FR][lev][2*n+1];
// 		HSBoxesLev[FR][lev][2*n]=BoxesLev[FR][lev][2*n]+zoom*Buffers[FR][2*n];
	      }
	  }
      }
    int RANKY;
    MPI_Comm_rank(FractalWorld,&RANKY);
    FRBoxesLev=BoxesLev[RANKY];
    FRBBoxesLev=BBoxesLev[RANKY];
    FRPBoxesLev=PBoxesLev[RANKY];
  }
  void Fractal_Memory::calc_RealBoxes()
  {
    //      cerr << "real " << FractalNodes << " " << grid_length << "\n";
    RealBoxes.resize(FractalNodes);
    RealPBoxes.resize(FractalNodes);
    RealIBoxes.resize(FractalNodes);
    LeftCorners.resize(FractalNodes);
    double glinv=1.0/static_cast<double>(grid_length);
    double DB=glinv;
    if(periodic)
      DB=-glinv;
    BigBox.resize(6);
    BigBox[0]=DB;
    BigBox[2]=DB;
    BigBox[4]=DB;
    BigBox[1]=1.0-DB;
    BigBox[3]=1.0-DB;
    BigBox[5]=1.0-DB;
    for(int b=0;b<FractalNodes;b++)
      {
	RealBoxes[b].resize(6);
	RealPBoxes[b].resize(6);
	RealIBoxes[b].resize(6);
	LeftCorners[b].resize(3);
	for(int ni=0;ni<6;ni+=2)
	  {
	    //	      cerr << " b ni " << b << " " << ni << "\n";
	    RealBoxes[b][ni]=static_cast<double>(Boxes[b][ni])*glinv;
	    RealBoxes[b][ni+1]=static_cast<double>(Boxes[b][ni+1]+1)*glinv;
	    RealPBoxes[b][ni]=static_cast<double>(PBoxes[b][ni])*glinv;
	    RealPBoxes[b][ni+1]=static_cast<double>(PBoxes[b][ni+1])*glinv;
	    LeftCorners[b][ni/2]=Boxes[b][ni];
	    RealIBoxes[b][ni]=RealBoxes[b][ni]+glinv*2.0;
	    RealIBoxes[b][ni+1]=RealBoxes[b][ni+1]-glinv*2.0;
	    if(periodic)
	      {
		if(Boxes[b][ni] == 0)
		  {
		    LeftCorners[b][ni/2]=-2;
		    RealIBoxes[b][ni]=-10.0;
		  }
		if(Boxes[b][ni+1]==grid_length-1)
		  RealIBoxes[b][ni+1]=10.0;
		continue;
	      }
	    else
	      {
		/*    I changed this */
		RealBoxes[b][ni]=max(RealBoxes[b][ni],glinv);
		RealBoxes[b][ni+1]=min(RealBoxes[b][ni+1],1.0-glinv);
		RealPBoxes[b][ni]=max(RealPBoxes[b][ni],glinv);
		RealPBoxes[b][ni+1]=min(RealPBoxes[b][ni+1],1.0-glinv);
	      }
	  }
      }
    //      cerr << " real b " << "\n";
  }
  int Fractal_Memory::fftw_where(int i,int j,int k,int lb,int lc) const
  {
    return k+(j+(i-p_mess->start_x)*lb)*lc;
  }
  void Fractal_Memory::Full_Stop() const
  {
    p_mess->Full_Stop();
  }
  void Fractal_Memory::make_scaling()
  {
    scaling=1.0;
    if(spectrum_number==1)
      {
	double a1=pow(46.9*omega_0*h*h,0.67)*(1.0+pow(32.1*omega_0*h*h,-0.532));
	double a2=pow(12.0*omega_0*h*h,0.424)*(1.0+pow(45.0*omega_0*h*h,-0.582));
	double alpha=pow(a1,-omega_b/omega_0)*pow(a2,-pow(omega_b/omega_0,3));
	scaling=box_length*omega_0*h*h*sqrt(alpha)*pow(1.0-omega_b/omega_0,0.6);
	cerr << "scaling " << a1 << " " << a2 << " " << alpha << " " << " " << box_length << " " << h << " " << scaling << "\n";
      }
  }
}
