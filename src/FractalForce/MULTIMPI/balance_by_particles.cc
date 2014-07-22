#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void balance_by_particles(Fractal_Memory* PFM)
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
    double scalepoint=1.0e-10;
    double scalepart=1.0;
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
	binary_balancing(numbersz,minimumz,FractalNodes2,real_length,lowerz,upperz);
      }
    PFM->p_mess->Send_INT_from_ROOT(lowerz,FractalNodes2,ROOTZ);
    PFM->p_mess->Send_INT_from_ROOT(upperz,FractalNodes2,ROOTZ);
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
	//	cout << " lower upper Z " << FractalRank << " " << FRZ << " " << lowerz[FRZ] << " " << upperz[FRZ] << "\n" ;
      }
    vector <int>numbert(FractalNodes2,0);
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
	numbert[FRZ]++;
      }
    /*
    cout << " numberT " << FractalRank << " ";
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      cout << FRZ << " " << numbert[FRZ];
    cout << "" << endl;
    */
    const int ROOTY=ROOTZ-FractalNodes2/2;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      PFM->p_mess->Find_Sum_DOUBLE_to_ROOT(numbersy[FRZ],real_length,ROOTY+FRZ);


    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	if(FractalRank == ROOTY+FRZ)
	  {
	    double minimumy=static_cast<double>(upperz[FRZ]-lowerz[FRZ])*real_length*scalepoint;
	    binary_balancing(numbersy[FRZ],minimumy,FractalNodes1,real_length,lowery[FRZ],uppery[FRZ]);
	  }
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	PFM->p_mess->Send_INT_from_ROOT(lowery[FRZ],FractalNodes1,ROOTY+FRZ);
	PFM->p_mess->Send_INT_from_ROOT(uppery[FRZ],FractalNodes1,ROOTY+FRZ);
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
	    //	    cout << " lower upper Y Z " << FractalRank << " " << FRY << " " << FRZ << " " << lowery[FRZ][FRY] << " " << uppery[FRZ][FRY];
	    //	    cout << " "   << lowerz[FRZ] << " " << upperz[FRZ] << "\n" ;
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
	    if(FractalRank == ROOTX+FRZ*FractalNodes1+FRY)
	      {
		double minimumx=static_cast<double>((upperz[FRZ]-lowerz[FRZ])*(uppery[FRZ][FRY]-lowery[FRZ][FRY]))*scalepoint;
		binary_balancing(numbersx[FRZ][FRY],minimumx,FractalNodes0,real_length,lowerx[FRZ][FRY],upperx[FRZ][FRY]);
	      }
	  }
      }
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    PFM->p_mess->Send_INT_from_ROOT(lowerx[FRZ][FRY],FractalNodes0,ROOTX+FRZ*FractalNodes1+FRY);
	    PFM->p_mess->Send_INT_from_ROOT(upperx[FRZ][FRY],FractalNodes0,ROOTX+FRZ*FractalNodes1+FRY);
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      alowerx[FRZ][FRY][FRX]=lowerx[FRZ][FRY][FRX];
	  }
      }

    double time4=PFM->p_mess->Clock();
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
		//		cout << " lower upper X Y Z " << FractalRank << " " << FRX << " " << FRY << " " << FRZ << " " << lowerx[FRZ][FRY][FRX] << " " << upperx[FRZ][FRY][FRX];
		//		cout << " " << lowery[FRZ][FRY] << " " << uppery[FRZ][FRY] << " " << lowerz[FRZ] << " " << upperz[FRZ] << "\n" ;
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
		fprintf(PFFM," BOXES %d \t %d \t %d \t %d \t %d ",FR,FRX,FRY,FRZ,VOL);
		fprintf(PFFM," \t %d \t %d \t %d ",PFM->Boxes[FR][0],PFM->Boxes[FR][1],PFM->Boxes[FR][2]);
		fprintf(PFFM," \t %d \t %d \t %d \n",PFM->Boxes[FR][3],PFM->Boxes[FR][4],PFM->Boxes[FR][5]);
		FR++;
	      }
	  }
      }
    double time5=PFM->p_mess->Clock();
    PFM->p_file->note(true," made new Boxes with equal smart particles ");
    PFM->calc_Buffers_and_more();
    PFM->p_file->note(true," made new Buffers with equal smart particles ");
    PFM->calc_RealBoxes();
    PFM->p_file->note(true," made new RealBoxes with equal smart particles ");
    PF->redo(PFM);
    PFM->p_file->note(true," redo fractal ");
    double time6=PFM->p_mess->Clock();
    fprintf(PFM->p_file->PFTime," balance particles %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E \n",time1-time0,time2-time1,time3-time2,time4-time3,time5-time4,time6-time5,time6-time1);
  }
}
