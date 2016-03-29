#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void balance_by_particles(Fractal_Memory* PFM)
  {
    double time0=PFM->p_mess->Clock();
    Fractal* PF=PFM->p_fractal;
    FILE* PFFM=PFM->p_file->PFFractalMemory;
    int FractalRank=PFM->p_mess->FractalRank;
    const int FractalNodes0=PFM->FractalNodes0;
    const int FractalNodes1=PFM->FractalNodes1;
    const int FractalNodes2=PFM->FractalNodes2;
    const int FractalNodes=PFM->FractalNodes;
    const double aFractalNodes0=FractalNodes0;
    const double aFractalNodes1=FractalNodes1;
    const double aFractalNodes2=FractalNodes2;
    //    int balance=PFM->balance;
    double scalepoint=0.0001;
    double scalepart=1.0;
    //    double scalepart=1.0e-6;
    int steps=PFM->steps;
    fprintf(PFFM," scalings balance %d %10.2E %10.2E \n",steps,scalepoint,scalepart);
    int real_length=PFM->grid_length;
    double alength=real_length;
    vector <double> numbersz(real_length,0.0);
    vector <double> pos(3);
    for(int ni=0;ni < PF->get_number_particles();++ni)
      {
	PF->particle_list[ni]->get_pos(pos);
	if(pos[0] < 0.0 || pos[0] >= 1.0 ||
	   pos[1] < 0.0 || pos[1] >= 1.0 ||
	   pos[2] < 0.0 || pos[2] >= 1.0)
	  continue;
	int nz=pos[2]*alength;
	assert(nz >= 0);
	assert(nz < real_length);
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
    const int ROOTZ=FractalNodes/2;
    PFM->p_mess->Find_Sum_DOUBLE_to_ROOT(numbersz,real_length,ROOTZ);
    if(FractalRank == ROOTZ)
      {
	double minimumz=alength*alength*scalepoint;
	vector <double>numbersz0(real_length);
	vector <double>snumbersz(real_length+1,0.0);
	numbersz0=numbersz;
	snumbersz[0]=0.0;
	for(int nz=0;nz<real_length;nz++)
	  snumbersz[nz+1]=snumbersz[nz]+numbersz0[nz]+minimumz; 
	
	double pztotal=snumbersz[real_length];
	for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
	  {
	    double aFRZ=FRZ;
	    double target=1.0e-30+(aFRZ*pztotal)/aFractalNodes2;
	    lowerz[FRZ]=std::lower_bound(snumbersz.begin(),snumbersz.end(),target)-snumbersz.begin()-1;
	    if(FRZ > 0 && lowerz[FRZ] <= lowerz[FRZ-1])
	      lowerz[FRZ]=lowerz[FRZ-1]+1;
	    fprintf(PFFM," target a %d \t %10.2E \t %d \t %d \n",FRZ,target,pztotal,lowerz[FRZ]);
	    if(FRZ > 0)
	      upperz[FRZ-1]=lowerz[FRZ];
	    alowerz[FRZ]=lowerz[FRZ];	  
	  }
	upperz[FractalNodes2-1]=real_length;
      }
    PFM->p_mess->Send_INT_from_ROOT(lowerz,FractalNodes2,ROOTZ);
    PFM->p_mess->Send_INT_from_ROOT(upperz,FractalNodes2,ROOTZ);
    double time1=PFM->p_mess->Clock();

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
	cout << " target b " << FractalRank << " " << FRZ << " " << lowerz[FRZ] << " " << upperz[FRZ] << "\n" ;
      }
    for(int ni=0;ni < PF->get_number_particles();++ni)
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
      }
    
    const int ROOTY=ROOTZ-FractalNodes2/2;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      PFM->p_mess->Find_Sum_DOUBLE_to_ROOT(numbersy[FRZ],real_length,ROOTY+FRZ);


    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	if(FractalRank != ROOTY+FRZ)
	  continue;
	double minimumy=static_cast<double>(upperz[FRZ]-lowerz[FRZ])*real_length*scalepoint;
	vector <double>numbersy0(real_length);
	vector <double>snumbersy(real_length+1,0.0);
	numbersy0=numbersy[FRZ];
	snumbersy[0]=0.0;
	for(int ny=0;ny<real_length;ny++)
	  snumbersy[ny+1]=snumbersy[ny]+numbersy0[ny]+minimumy; 
	
	double pytotal=snumbersy[real_length];
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    double aFRY=FRY;
	    double target=1.0e-30+(aFRY*pytotal)/aFractalNodes1;
	    lowery[FRZ][FRY]=std::lower_bound(snumbersy.begin(),snumbersy.end(),target)-snumbersy.begin()-1;
	    if(FRY > 0 && lowery[FRZ][FRY] <= lowery[FRZ][FRY-1])
	      lowery[FRZ][FRY]=lowery[FRZ][FRY-1]+1;
	    fprintf(PFFM," target b %d \t %10.2E \t %d \t %d \n",FRY,target,pytotal,lowery[FRZ][FRY]);
	    if(FRY > 0)
	      uppery[FRZ][FRY-1]=lowery[FRZ][FRY];
	    alowery[FRZ][FRY]=lowery[FRZ][FRY];
	  }
	uppery[FRZ][FractalNodes1-1]=real_length;
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	PFM->p_mess->Send_INT_from_ROOT(lowery[FRZ],FractalNodes1,ROOTY+FRZ);
	PFM->p_mess->Send_INT_from_ROOT(uppery[FRZ],FractalNodes1,ROOTY+FRZ);
      }

    double time2=PFM->p_mess->Clock();
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    lowerx[FRZ][FRY].resize(FractalNodes0);
	    alowerx[FRZ][FRY].resize(FractalNodes0);
	    upperx[FRZ][FRY].resize(FractalNodes0);
	    numbersx[FRZ][FRY].assign(real_length,0.0);
	    cout << " target c " << FractalRank << " " << FRY << " " << FRZ << " " << lowery[FRZ][FRY] << " " << lowerz[FRZ];
	    cout << " "   << uppery[FRZ][FRY] << " " << upperz[FRZ] << "\n" ;
	  }
      }
    for(int ni=0;ni < PF->get_number_particles();++ni)
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
    const int ROOTX=ROOTZ-FractalNodes2*FractalNodes1/2;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      for(int FRY=0;FRY<FractalNodes1;FRY++)
	PFM->p_mess->Find_Sum_DOUBLE_to_ROOT(numbersx[FRZ][FRY],real_length,ROOTX+FRZ*FractalNodes1+FRY);

    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    if(FractalRank != ROOTX+FRZ*FractalNodes1+FRY)
	      continue;
	    double minimumx=static_cast<double>((upperz[FRZ]-lowerz[FRZ])*(uppery[FRZ][FRY]-lowery[FRZ][FRY]))*scalepoint;
	    vector <double>numbersx0(real_length);
	    vector <double>snumbersx(real_length+1,0.0);
	    numbersx0=numbersx[FRZ][FRY];
	    snumbersx[0]=0.0;
	    for(int nx=0;nx<real_length;nx++)
	      snumbersx[nx+1]=snumbersx[nx]+numbersx0[nx]+minimumx; 
	    
	    double pxtotal=snumbersx[real_length];
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		double aFRX=FRX;
		double target=1.0e-30+(aFRX*pxtotal)/aFractalNodes0;
		lowerx[FRZ][FRY][FRX]=std::lower_bound(snumbersx.begin(),snumbersx.end(),target)-snumbersx.begin()-1;
		if(FRX > 0 && lowerx[FRZ][FRY][FRX] <= lowerx[FRZ][FRY][FRX-1])
		  lowerx[FRZ][FRY][FRX]=lowerx[FRZ][FRY][FRX-1]+1;
		fprintf(PFFM," target %d \t %10.2E \t %d \t %d \n",FRX,target,pxtotal,lowerx[FRZ][FRY][FRX]);
		if(FRX > 0)
		  upperx[FRZ][FRY][FRX-1]=lowerx[FRZ][FRY][FRX];
		alowerx[FRZ][FRY][FRX]=lowerx[FRZ][FRY][FRX];
	      }
	    upperx[FRZ][FRY][FractalNodes0-1]=real_length;
	  }
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    PFM->p_mess->Send_INT_from_ROOT(lowerx[FRZ][FRY],FractalNodes0,ROOTX+FRZ*FractalNodes1+FRY);
	    PFM->p_mess->Send_INT_from_ROOT(upperx[FRZ][FRY],FractalNodes0,ROOTX+FRZ*FractalNodes1+FRY);
	  }
      }

    double time3=PFM->p_mess->Clock();
    int FR=0;
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
		cout << " target d " << FractalRank << " " << FRX << " " << FRY << " " << FRZ << " " << lowerx[FRZ][FRY][FRX] << " " << lowery[FRZ][FRY];
		cout << " " << lowerz[FRZ] << " " << upperx[FRZ][FRY][FRX] << " " << uppery[FRZ][FRY] << " " << upperz[FRZ] << "\n" ;
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
		VOL*=(PFM->Boxes[FR][3]-PFM->Boxes[FR][2]);
		VOL*=(PFM->Boxes[FR][5]-PFM->Boxes[FR][4]);
		fprintf(PFFM," BOXES %d \t %d \t %d \t %d \t %d ",FR,FRX,FRY,FRZ,VOL);
		fprintf(PFFM," %d \t %d \t %d ",PFM->Boxes[FR][0],PFM->Boxes[FR][1],PFM->Boxes[FR][2]);
		fprintf(PFFM," \t %d \t %d \t %d \n",PFM->Boxes[FR][3],PFM->Boxes[FR][4],PFM->Boxes[FR][5]);
		FR++;
	      }
	  }
      }
    double time4=PFM->p_mess->Clock();
    PFM->p_file->note(true," made new Boxes with equal smart particles ");
    PFM->calc_Buffers_and_more();
    PFM->p_file->note(true," made new Buffers with equal smart particles ");
    PFM->calc_RealBoxes();
    PFM->p_file->note(true," made new RealBoxes with equal smart particles ");
    PF->redo(PFM);
    PFM->p_file->note(true," redo fractal ");
    double time5=PFM->p_mess->Clock();
    fprintf(PFM->p_file->PFTime," balance particles %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E \n",time1-time0,time2-time1,time3-time2,time4-time3,time5-time4,time5-time0);
  }
}
