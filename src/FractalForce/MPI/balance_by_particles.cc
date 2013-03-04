#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void balance_by_particles(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    ofstream& FFM=PFM->p_file->FileFractalMemory;
    const int ROOT=0;
    const int mult=10;

    const int FractalNodes0=PFM->FractalNodes0;
    const int FractalNodes1=PFM->FractalNodes1;
    const int FractalNodes2=PFM->FractalNodes2;

    const int length0=FractalNodes0*mult;
    const int length1=FractalNodes1*mult;
    const int length2=FractalNodes2*mult;
    const int length3=length0*length1*length2;

    const double alength0=length0;
    const double alength1=length1;
    const double alength2=length2;
    const double alength3=length3;

    const int real_length=PFM->grid_length;
    const int real_cells=Misc::pow(real_length,3);
    const double a_real_cells=real_cells;
    const int min_count=a_real_cells/alength3;
    int* numbers=new int[length3];
    for(int ni=0;ni<length3;ni++)
      numbers[ni]=0;
    vector <double> pos(3);
    for(int ni=0;ni < PF->get_number_particles();++ni)
      {
	PF->particle_list[ni]->get_pos(pos);
	if(pos[0] < 0.0 || pos[0] >= 1.0 ||
	   pos[1] < 0.0 || pos[1] >= 1.0 ||
	   pos[2] < 0.0 || pos[2] >= 1.0)
	  continue;
	int nx=pos[0]*alength0;
	int ny=pos[1]*alength1;
	int nz=pos[2]*alength2;
	assert(nx >= 0 && nx < alength0);
	assert(ny >= 0 && ny < alength1);
	assert(nz >= 0 && nz < alength2);
	int n=nx+(ny+nz*length1)*length0;
	numbers[n]++;
      }
    PFM->p_mess->Find_Sum_INT_to_ROOT(numbers,length3,ROOT);
    PFM->p_mess->Send_INT_from_ROOT(numbers,length3,ROOT);
    vector <int>pointsa(length3);
    for(int ni=0;ni<length3;ni++)
      pointsa[ni]=numbers[ni]+min_count;
    delete [] numbers;
    numbers=0;
    vector <int>sumz(length2+1,0);
    for(int nz=0;nz<length2;nz++)
      {
	Misc::sum_up(sumz[nz+1],pointsa,nz*length0*length1,(nz+1)*length0*length1,1);
	sumz[nz+1]+=sumz[nz];
	FFM << " sumz " << nz << " " << sumz[nz+1] << endl;
      }
    vector <int>lowerz(FractalNodes2);
    vector <int>upperz(FractalNodes2);
    vector < vector <int> > lowery(FractalNodes2);
    vector < vector <int> > uppery(FractalNodes2);
    vector < vector < vector <int> > > lowerx(FractalNodes2);
    vector < vector < vector <int> > > upperx(FractalNodes2);
    int pztotal=sumz[length2];
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	int target=(FRZ*pztotal)/FractalNodes2;
	lowerz[FRZ]=std::lower_bound(sumz.begin(),sumz.end(),target)-sumz.begin();
	FFM << " target " << FRZ << " " << target << " " << pztotal << " " << lowerz[FRZ] << endl;
	if(FRZ > 0)
	  upperz[FRZ-1]=lowerz[FRZ];
      }
    upperz[FractalNodes2-1]=length2;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	FFM << " FRZ " << FRZ << " " << lowerz[FRZ] << " " << upperz[FRZ] << " ";
	FFM << sumz[lowerz[FRZ]] << " " << sumz[upperz[FRZ]] << endl;
      }

    vector <int>pointsb(length0*length1,0);
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	lowery[FRZ].resize(FractalNodes1);
	uppery[FRZ].resize(FractalNodes1);
	lowerx[FRZ].resize(FractalNodes1);
	upperx[FRZ].resize(FractalNodes1);
	int ni=0;
	int stride=length0*length1;
	int start=stride*lowerz[FRZ];
	int end=stride*upperz[FRZ];
	for(int ny=0;ny<length1;ny++)
	  {
	    for(int nx=0;nx<length0;nx++)
	      {
		Misc::sum_up(pointsb[ni],pointsa,ni+start,ni+end,stride);
		ni++;
	      }
	  }
	vector <int>sumy(length1+1,0);
	for(int ny=0;ny<length1;ny++)
	  {
	    Misc::sum_up(sumy[ny+1],pointsb,ny*length0,(ny+1)*length0,1);
	    sumy[ny+1]+=sumy[ny];
	  }
	int pytotal=sumy[length1];
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    int target=(FRY*pytotal)/FractalNodes1;
	    lowery[FRZ][FRY]=std::lower_bound(sumy.begin(),sumy.end(),target)-sumy.begin();
	    if(FRY > 0)
	      uppery[FRZ][FRY-1]=lowery[FRZ][FRY];
	  }
	uppery[FRZ][FractalNodes1-1]=length1;
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    FFM << " FRYZ " << FRZ << " " << FRY << " " << lowery[FRZ][FRY] << " " << uppery[FRZ][FRY] << " ";
	    FFM << sumy[lowery[FRZ][FRY]] << " " << sumy[uppery[FRZ][FRY]] << endl;
	  }

	vector <int>pointsc(length0,0);
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    lowerx[FRZ][FRY].resize(FractalNodes0);
	    upperx[FRZ][FRY].resize(FractalNodes0);
	    int stride=length0;
	    int start=lowery[FRZ][FRY]*length0;
	    int end=uppery[FRZ][FRY]*length0;
	    for(int nx=0;nx<length0;nx++)
	      {
		Misc::sum_up(pointsc[nx],pointsb,nx+start,nx+end,stride);
	      }
	    vector <int>sumx(length0+1,0);
	    for(int nx=0;nx<length0;nx++)
	      sumx[nx+1]=sumx[nx]+pointsc[nx];
	    int pxtotal=sumx[length0];
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		int target=(FRX*pxtotal)/FractalNodes0;
		lowerx[FRZ][FRY][FRX]=std::lower_bound(sumx.begin(),sumx.end(),target)-sumx.begin();
		if(FRX > 0)
		  upperx[FRZ][FRY][FRX-1]=lowerx[FRZ][FRY][FRX];
	      }
	    upperx[FRZ][FRY][FractalNodes0-1]=length0;
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		FFM << " FRXYZ " << FRZ << " " << FRY << " " << FRX << " " << lowerx[FRZ][FRY][FRX] << " " << upperx[FRZ][FRY][FRX] << " " ;
		FFM << sumx[lowerx[FRZ][FRY][FRX]] << " " << sumx[upperx[FRZ][FRY][FRX]] << endl;
	      }
	  }
      }
    //    PFM->p_mess->Full_Stop();
    //    assert(0);
    
    int FR=0;
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	int BZ0=(lowerz[FRZ]*real_length)/length2;
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    int BY0=(lowery[FRZ][FRY]*real_length)/length1;
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		int BX0=(lowerx[FRZ][FRY][FRX]*real_length)/length0;
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
		FFM << " BOXES " << FR << " " << PFM->Boxes[FR][0] << " " << PFM->Boxes[FR][1] << " " << PFM->Boxes[FR][2];
		FFM << " " << PFM->Boxes[FR][3] << " " << PFM->Boxes[FR][4] << " " << PFM->Boxes[FR][5] << endl;
		FR++;
	      }
	  }
      }
    PFM->p_file->note(true," made new Boxes with equal particles ");
    PFM->calc_Buffers_and_more();
    PFM->p_file->note(true," made new Buffers with equal particles ");
    PFM->calc_RealBoxes();
    PFM->p_file->note(true," made new RealBoxes with equal particles ");
    PF->redo(PFM);
    PFM->p_file->note(true," redo fractal ");
  }
}
