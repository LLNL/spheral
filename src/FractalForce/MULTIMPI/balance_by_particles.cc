#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void balance_by_particles(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    ofstream& FFM=PFM->p_file->FileFractalMemory;
    int balance=PFM->balance;
    double scalepoint=1.0;
    double scalepart=7.0;
    int steps=PFM->steps;
    /*
    if(steps % 20 < 10 )
      {
	scalepoint= 1.0+1.5*static_cast<double>(steps % 10);
	scalepart= 1.0;
      }
    else
      {
	scalepoint=1.0;
	scalepart= 1.0+1.5*static_cast<double>(steps % 10);
      }
    */
    FFM << " scalings balance " << steps << " " << scalepoint << " " << scalepart << endl;
    const int ROOT=0;
    const int mult=10;
    const double MN=PFM->minimum_number;
    const double logMNinv=1.0/log(MN);

    const int FractalNodes0=PFM->FractalNodes0;
    const int FractalNodes1=PFM->FractalNodes1;
    const int FractalNodes2=PFM->FractalNodes2;
    const double aFractalNodes0=FractalNodes0;
    const double aFractalNodes1=FractalNodes1;
    const double aFractalNodes2=FractalNodes2;

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
    const double min_count=a_real_cells/alength3;
    const double min_count_inv=1.0/min_count;
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
	assert(nx >= 0 && nx < length0);
	assert(ny >= 0 && ny < length1);
	assert(nz >= 0 && nz < length2);
	int n=nx+(ny+nz*length1)*length0;
	numbers[n]++;
      }
    PFM->p_mess->Find_Sum_INT_to_ROOT(numbers,length3,ROOT);
    PFM->p_mess->Send_INT_from_ROOT(numbers,length3,ROOT);
    vector <double>pointsa(length3);
    for(int ni=0;ni<length3;ni++)
      pointsa[ni]=numbers[ni];
    delete [] numbers;
    numbers=0;
    for(int ni=0;ni<length3;ni++)
      {
	if(balance > 1)
	  {
	    double lev=log(pointsa[ni]*min_count_inv)*logMNinv;
	    pointsa[ni]*=(lev+1.0)*scalepart;
	    pointsa[ni]+=min_count*pow(8.0,lev)*scalepoint;
	  }
	else
	  {
	    pointsa[ni]=pointsa[ni]*scalepart+min_count*scalepoint;
	  }

      }
    vector <double>sumz(length2+1,0);
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
    double pztotal=sumz[length2];
    for(int FRZ=0;FRZ<FractalNodes2;FRZ++)
      {
	double aFRZ=FRZ;
	double target=(aFRZ*pztotal)/aFractalNodes2;
	lowerz[FRZ]=std::lower_bound(sumz.begin(),sumz.end(),target)-sumz.begin();
	if(FRZ > 0 && lowerz[FRZ] <= lowerz[FRZ-1])
	  lowerz[FRZ]=lowerz[FRZ-1]+1;
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

    vector <double>pointsb(length0*length1,0);
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
	vector <double>sumy(length1+1,0);
	for(int ny=0;ny<length1;ny++)
	  {
	    Misc::sum_up(sumy[ny+1],pointsb,ny*length0,(ny+1)*length0,1);
	    sumy[ny+1]+=sumy[ny];
	  }
	double pytotal=sumy[length1];
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    double aFRY=FRY;
	    double target=(aFRY*pytotal)/aFractalNodes1;
	    lowery[FRZ][FRY]=std::lower_bound(sumy.begin(),sumy.end(),target)-sumy.begin();
	    if(FRY > 0 && lowery[FRZ][FRY] <= lowery[FRZ][FRY-1])
	      lowery[FRZ][FRY]=lowery[FRZ][FRY-1]+1;
	    if(FRY > 0)
	      uppery[FRZ][FRY-1]=lowery[FRZ][FRY];
	  }
	uppery[FRZ][FractalNodes1-1]=length1;
	for(int FRY=0;FRY<FractalNodes1;FRY++)
	  {
	    FFM << " FRYZ " << FRZ << " " << FRY << " " << lowery[FRZ][FRY] << " " << uppery[FRZ][FRY] << " ";
	    FFM << sumy[lowery[FRZ][FRY]] << " " << sumy[uppery[FRZ][FRY]] << endl;
	  }

	vector <double>pointsc(length0,0);
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
	    vector <double>sumx(length0+1,0);
	    for(int nx=0;nx<length0;nx++)
	      sumx[nx+1]=sumx[nx]+pointsc[nx];
	    double pxtotal=sumx[length0];
	    for(int FRX=0;FRX<FractalNodes0;FRX++)
	      {
		double aFRX=FRX;
		double target=(aFRX*pxtotal)/aFractalNodes0;
		lowerx[FRZ][FRY][FRX]=std::lower_bound(sumx.begin(),sumx.end(),target)-sumx.begin();
		if(FRX > 0 && lowerx[FRZ][FRY][FRX] <= lowerx[FRZ][FRY][FRX-1])
		  lowerx[FRZ][FRY][FRX]=lowerx[FRZ][FRY][FRX-1]+1;
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
		int VOL=(PFM->Boxes[FR][1]-PFM->Boxes[FR][0]);
		VOL*=(PFM->Boxes[FR][3]-PFM->Boxes[FR][2]);
		VOL*=(PFM->Boxes[FR][5]-PFM->Boxes[FR][4]);
		FFM << " BOXES " << FR << " " << FRX << " " << FRY << " " << FRZ << " " << VOL << " ";
		FFM << PFM->Boxes[FR][0] << " " << PFM->Boxes[FR][1] << " " << PFM->Boxes[FR][2];
		FFM << " " << PFM->Boxes[FR][3] << " " << PFM->Boxes[FR][4] << " " << PFM->Boxes[FR][5] << endl;
		FR++;
	      }
	  }
      }
    PFM->p_file->note(true," made new Boxes with equal smart particles ");
    PFM->calc_Buffers_and_more();
    PFM->p_file->note(true," made new Buffers with equal smart particles ");
    PFM->calc_RealBoxes();
    PFM->p_file->note(true," made new RealBoxes with equal smart particles ");
    PF->redo(PFM);
    PFM->p_file->note(true," redo fractal ");
  }
}
