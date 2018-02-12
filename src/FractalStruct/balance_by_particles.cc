#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void balance_by_particles(Fractal_Memory* PFM,bool withparts)
  {
    double time0=PFM->p_mess->Clock();
    PFM->p_mess->Full_Stop();
    double time1=PFM->p_mess->Clock();
    Fractal* PF=PFM->p_fractal;
    FILE* PFFM=PFM->p_file->PFFractalMemory;
    int FractalRank=PFM->p_mess->FractalRank;
    const int FractalNodes0=PFM->FractalNodes0;
    const int FractalNodes1=PFM->FractalNodes1;
    const int FractalNodes2=PFM->FractalNodes2;
    const int FractalNodes=PFM->FractalNodes;
    const int FFTNodes=PFM->FFTNodes;
    const int ROOTNODE=PFM->p_mess->ROOTNODE;
    int how_many_particles=0;
    if(withparts)
      how_many_particles=PF->get_number_particles();
    double NOSHRINK=1.0;
    double SHRINK=0.2;
    if(PFM->p_mess->WallNodes < FFTNodes)
      SHRINK=NOSHRINK;

    SHRINK=NOSHRINK;     //Shrinking of FFT nodes turned off for now.

    vector <double>targets(FractalNodes2+1);
    targets[0]=0.0;
    int FR=0;
    for(int FR2=0;FR2<FractalNodes2;FR2++)
      {
	targets[FR2+1]=targets[FR2];
	for(int FR1=0;FR1<FractalNodes1;FR1++)
	  {
	    for(int FR0=0;FR0<FractalNodes0;FR0++)
	      {
		if(PFM->p_mess->ItIsAnFFTNode[FR])
		  targets[FR2+1]+=SHRINK;
		else
		  targets[FR2+1]+=NOSHRINK;
		FR++;
	      }
	  }
      }
    for(int FR2=0;FR2<=FractalNodes2;FR2++)
      targets[FR2]/=targets[FractalNodes2];
    double scalepoint=1.0e-10;
    double scalepart=1.0;
    int steps=PFM->steps;
    if(FractalRank == 0)
      fprintf(PFFM," scalings balance %d %10.2E %10.2E \n",steps,scalepoint,scalepart);
    int real_length=PFM->grid_length;
    double alength=real_length;
    vector <double> numbersz(real_length,0.0);
    vector <double> pos(3);
    bool OKA=true;
    bool OKB=true;
    for(int ni=0;ni < how_many_particles;++ni)
      {
	PF->particle_list[ni]->get_pos(pos);
	if(pos[0] < 0.0 || pos[0] >= 1.0 ||
	   pos[1] < 0.0 || pos[1] >= 1.0 ||
	   pos[2] < 0.0 || pos[2] >= 1.0)
	  continue;
	int nz=pos[2]*alength;
	if(nz < 0 || nz >= real_length)
	  {
	    fprintf(PFFM,"ERROR in balanceA %13.5E %13.5E %13.5E %d %d %d %d \n",pos[0],pos[1],pos[2],alength,ni,nz,real_length);
	    cerr << "ERROR in balanceA %13.5E %13.5E %13.5E %d %d %d %d "<< FractalRank << " " << pos[0]<< " " << pos[1]<< " " << pos[2]<< " " << alength<< " " << ni<< " " << nz<< " " << real_length << endl;
	    if(nz < 0)
	      OKA=false;;
	    if(nz >= real_length)
	      OKB=false;
	    nz=max(nz,0);
	    nz=min(nz,real_length-1);
	  }
	numbersz[nz]+=scalepart;
      }

    vector <int>lowerz(FractalNodes2);
    vector <double>alowerz(FractalNodes2);
    vector <int>upperz(FractalNodes2);
    vector < vector <int> > lowery(FractalNodes2);
    vector < vector <double> > alowery(FractalNodes2);
    vector < vector <int> > uppery(FractalNodes2);
    vector < vector < vector <int> > > lowerx(FractalNodes2);
    vector < vector < vector <double> > > alowerx(FractalNodes2);
    vector < vector < vector <int> > > upperx(FractalNodes2);
    vector < vector <double> > numbersy(FractalNodes2);
    vector < vector < vector <double> > > numbersx(FractalNodes2);
    PFM->p_mess->Find_Sum_DOUBLE_to_ROOT(numbersz,real_length,ROOTNODE);
    assert(OKA);
    assert(OKB);
    if(FractalRank == ROOTNODE)
      {
	double minimumz=alength*alength*scalepoint;
	binary_balancing(PFM,numbersz,minimumz,FractalNodes2,real_length,targets,lowerz,upperz);
      }
    PFM->p_mess->Send_INT_from_ROOT(lowerz,FractalNodes2,ROOTNODE);
    for(int FRZ=1;FRZ<FractalNodes2;FRZ++)
      upperz[FRZ-1]=lowerz[FRZ];
    upperz[FractalNodes2-1]=real_length-1;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      alowerz[FRZ]=lowerz[FRZ];	  
    double time2=PFM->p_mess->Clock();

    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	lowery[FRZ].resize(FractalNodes1);
	alowery[FRZ].resize(FractalNodes1);
	uppery[FRZ].resize(FractalNodes1);
	lowerx[FRZ].resize(FractalNodes1);
	alowerx[FRZ].resize(FractalNodes1);
	upperx[FRZ].resize(FractalNodes1);
	numbersy[FRZ].assign(real_length,0.0);
	numbersx[FRZ].resize(FractalNodes1);
      }
    vector <int>numbert(FractalNodes2,0);
    for(int ni=0;ni < how_many_particles;++ni)
      {
	PF->particle_list[ni]->get_pos(pos);
	if(pos[0] < 0.0 || pos[0] >= 1.0 ||
	   pos[1] < 0.0 || pos[1] >= 1.0 ||
	   pos[2] < 0.0 || pos[2] >= 1.0)
	  continue;
	double anz=1.0e-30+pos[2]*alength;
	int FRZ=std::lower_bound(alowerz.begin(),alowerz.end(),anz)-alowerz.begin()-1;
	assert(FRZ >= 0);
	assert(FRZ < FractalNodes2);
	int ny=pos[1]*alength;
	assert(ny >= 0);
	assert(ny < real_length);
	numbersy[FRZ][ny]+=scalepart;
	numbert[FRZ]++;
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	vector <float>numf(real_length);
	for(int ni=0;ni<real_length;ni++)
	  numf[ni]=numbersy[FRZ][ni];
	PFM->p_mess->Find_Sum_FLOAT_to_ROOT(numf,real_length,ROOTNODE);
	for(int ni=0;ni<real_length;ni++)
	  numbersy[FRZ][ni]=numf[ni];
      }
    targets.resize(FractalNodes1+1);

    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	int FR=FractalNodes0*FractalNodes1*FRZ;
	targets[0]=0.0;
	for(int FR1=0;FR1<FractalNodes1;FR1++)
	  {
	    targets[FR1+1]=targets[FR1];
	    for(int FR0=0;FR0<FractalNodes0;FR0++)
	      {
		if(PFM->p_mess->ItIsAnFFTNode[FR])
		  targets[FR1+1]+=SHRINK;
		else
		  targets[FR1+1]+=NOSHRINK;
		FR++;
	      }
	  }
	for(int FR1=0;FR1<=FractalNodes1;FR1++)
	  targets[FR1]/=targets[FractalNodes1];
	if(FractalRank == ROOTNODE)
	  {
	    double minimumy=static_cast<double>(upperz[FRZ]-lowerz[FRZ])*real_length*scalepoint;
	    binary_balancing(PFM,numbersy[FRZ],minimumy,FractalNodes1,real_length,targets,lowery[FRZ],uppery[FRZ]);
	  }
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	PFM->p_mess->Send_INT_from_ROOT(lowery[FRZ],FractalNodes1,ROOTNODE);
	for(int FRY=1;FRY<FractalNodes1;FRY++)
	  uppery[FRZ][FRY-1]=lowery[FRZ][FRY];
	uppery[FRZ][FractalNodes1-1]=real_length-1;
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  alowery[FRZ][FRY]=lowery[FRZ][FRY];
      }

    double time3=PFM->p_mess->Clock();
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    lowerx[FRZ][FRY].resize(FractalNodes0);
	    alowerx[FRZ][FRY].resize(FractalNodes0);
	    upperx[FRZ][FRY].resize(FractalNodes0);
	    numbersx[FRZ][FRY].assign(real_length,0.0);
	  }
      }
    for(int ni=0;ni < how_many_particles;++ni)
      {
	PF->particle_list[ni]->get_pos(pos);
	if(pos[0] < 0.0 || pos[0] >= 1.0 ||
	   pos[1] < 0.0 || pos[1] >= 1.0 ||
	   pos[2] < 0.0 || pos[2] >= 1.0)
	  continue;
	double anz=1.0e-30+pos[2]*alength;
	int FRZ=std::lower_bound(alowerz.begin(),alowerz.end(),anz)-alowerz.begin()-1;
	assert(FRZ >= 0);
	assert(FRZ < FractalNodes2);
	double any=1.0e-30+pos[1]*alength;
	int FRY=std::lower_bound(alowery[FRZ].begin(),alowery[FRZ].end(),any)-alowery[FRZ].begin()-1;
	assert(FRY >= 0);
	assert(FRY < FractalNodes1);
	int nx=pos[0]*alength;
	assert(nx >= 0);
	assert(nx < real_length);
	numbersx[FRZ][FRY][nx]+=scalepart;
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	vector <float>manynumbers;
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    for(int nx=0;nx<real_length;nx++)
	      manynumbers.push_back(numbersx[FRZ][FRY][nx]);
	  }
	PFM->p_mess->Find_Sum_FLOAT_to_ROOT(manynumbers,real_length*FractalNodes1,ROOTNODE);
	if(FractalRank == ROOTNODE)
	  {
	    int ni=0;
	    for(int FRY=0;FRY<FractalNodes1;FRY++)
	      {
		for(int nx=0;nx<real_length;nx++)
		  {
		    numbersx[FRZ][FRY][nx]=manynumbers[ni];
		    ni++;
		  }
	      }
	  }
      }
    targets.resize(FractalNodes0+1);
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    int FR=FRZ*FractalNodes0*FractalNodes1+FRY*FractalNodes0;
	    targets[0]=0.0;
	    for(int FR0=0;FR0<FractalNodes0;FR0++)
	      {
		if(PFM->p_mess->ItIsAnFFTNode[FR])
		  targets[FR0+1]=targets[FR0]+SHRINK;
		else
		  targets[FR0+1]=targets[FR0]+NOSHRINK;
		FR++;
	      }
	    for(int FR0=0;FR0<FractalNodes0;FR0++)
	      {
		targets[FR0]/=targets[FractalNodes0];
		//		if(FractalRank == 0)
		//		  cerr << " Target0 " << FR0 << " " << targets[FR0] << "\n";
	      }
	    if(FractalRank == ROOTNODE)
	      {
		double minimumx=static_cast<double>((upperz[FRZ]-lowerz[FRZ])*(uppery[FRZ][FRY]-lowery[FRZ][FRY]))*scalepoint;
		binary_balancing(PFM,numbersx[FRZ][FRY],minimumx,FractalNodes0,real_length,targets,lowerx[FRZ][FRY],upperx[FRZ][FRY]);
	      }
	  }
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    PFM->p_mess->Send_INT_from_ROOT(lowerx[FRZ][FRY],FractalNodes0,ROOTNODE);
	    for(int FRX=1;FRX<FractalNodes0;FRX++)
	      upperx[FRZ][FRY][FRX-1]=lowerx[FRZ][FRY][FRX];
	    upperx[FRZ][FRY][FractalNodes0-1]=real_length-1;
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      alowerx[FRZ][FRY][FRX]=lowerx[FRZ][FRY][FRX];
	  }
      }

    double time4=PFM->p_mess->Clock();
    FR=0;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	int BZ0=lowerz[FRZ];
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    int BY0=lowery[FRZ][FRY];
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		int BX0=lowerx[FRZ][FRY][FRX];
		PFM->Boxes[FR][0]=BX0;
		PFM->Boxes[FR][2]=BY0;
		PFM->Boxes[FR][4]=BZ0;
		FR++;
	      } 
	  }
      }
    FR=0;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		PFM->Boxes[FR][1]=real_length-1;
		if(FRX < FractalNodes0-1)
		  PFM->Boxes[FR][1]=PFM->Boxes[FR+1][0]-1;
		PFM->Boxes[FR][3]=real_length-1;
		if(FRY < FractalNodes1-1)
		  PFM->Boxes[FR][3]=PFM->Boxes[FR+FractalNodes0][2]-1;
		PFM->Boxes[FR][5]=real_length-1;
		if(FRZ < FractalNodes2-1)
		  PFM->Boxes[FR][5]=PFM->Boxes[FR+FractalNodes0*FractalNodes1][4]-1;
		int VOL=(PFM->Boxes[FR][1]-PFM->Boxes[FR][0]);
		assert(VOL > 0);
		VOL*=(PFM->Boxes[FR][3]-PFM->Boxes[FR][2]);
		assert(VOL > 0);
		VOL*=(PFM->Boxes[FR][5]-PFM->Boxes[FR][4]);
		assert(VOL > 0);
		if(FractalRank == 0)
		  {
		    fprintf(PFFM," BOXES %5d %5d %5d %5d %7d ",FR,FRX,FRY,FRZ,VOL);
		    fprintf(PFFM," %5d %5d %5d ",PFM->Boxes[FR][0],PFM->Boxes[FR][1],PFM->Boxes[FR][2]);
		    fprintf(PFFM," %5d %5d %5d \n",PFM->Boxes[FR][3],PFM->Boxes[FR][4],PFM->Boxes[FR][5]);
		  }
		FR++;
	      }
	  }
      }
    double time5=PFM->p_mess->Clock();
    if(FractalRank == 0)
      PFM->p_file->note(true," made new Boxes with equal smart particles ");
    PFM->calc_Buffers_and_more();
    if(FractalRank == 0)
      PFM->p_file->note(true," made new Buffers with equal smart particles ");
    PFM->calc_RealBoxes();
    if(FractalRank == 0)
      PFM->p_file->note(true," made new RealBoxes with equal smart particles ");
    PF->redo(PFM);
    if(FractalRank == 0)
      PFM->p_file->note(true," redo fractal ");
    double time6=PFM->p_mess->Clock();
    if(FractalRank == 0)
      fprintf(PFM->p_file->PFTime," balance particles %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E \n",time1-time0,time2-time1,time3-time2,time4-time3,time5-time4,time6-time5,time6-time1);
  }
}
