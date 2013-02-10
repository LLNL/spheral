#ifndef _Mess_Defined_
#define _Mess_Defined_
namespace FractalSpace
{
  class Mess{
  public:
    int FractalRank;
    int FractalNodes;
    int number_particles_total;
    Particle* parts_tmp;
    Particle* parts_interface;
    vector < vector <int> > Slices;
    vector < vector <int> > BoxS;
    vector < vector <int> > BoxSL;
    ptrdiff_t start_x;
    ptrdiff_t length_x;
    ptrdiff_t total_memory;
    fftw_plan plan_rc;
    fftw_plan plan_cr;
    double* potR;
    fftw_complex* potC; 
    vector <int> WhichSlice;
    vector <int>return_Slice_pos;
    vector <int>return_group;
    vector <int>return_point;
    vector <int>return_node;
    vector <int>what_Slice_point;
    vector <double> green;
    MPI_Comm FractalWorld;
    MPI_Comm HypreWorld;
    MPI_Comm SELF;
    Mess():
      FractalRank(0),
      FractalNodes(1),
      number_particles_total(-1),
      start_x(0),
      length_x(-1),
      total_memory(-1)
    {
      cout << " Empty Mess " << endl;
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP):
      FractalRank(0),
      FractalNodes(1),
      number_particles_total(-1),
      start_x(0),
      length_x(GR),
      total_memory(-1)
    {
      int grid_length=GR;
      bool periodic=PR;
      if(MR)
	{
	  MPIStartup();
	  FractalRank=what_is_my_rank(); 
	  FractalNodes=how_many_nodes(); 
	  FFTWStartup(grid_length,periodic);
	  calc_fftw_Slices(grid_length,periodic);	
	  calc_total_particles(NP);
	}
      else
	{
	  number_particles_total=NP;
	  length_x=grid_length;
	}
      cout << " made a mess " << FractalRank << " " << FractalNodes << " " << length_x << " " << start_x << " " << total_memory << endl;
    }
    ~Mess()
    {
      cout << " cleaned up a mess " << FractalRank << endl;
      FFTWFinal();
    }
    void MPIStartup()
    {
      int knights;
      MPI_Initialized(&knights);
      if(!knights)
	MPI_Init(NULL,NULL);
      FractalWorld=MPI_COMM_WORLD;
      HypreWorld=MPI_COMM_WORLD;
      SELF=MPI_COMM_SELF;
      FractalRank=what_is_my_rank(); 
      FractalNodes=how_many_nodes(); 
      cout << " initialized MPI " << FractalRank << " " << FractalNodes << endl;
    }
    void MPIFinal()
    {
      int knights;
      MPI_Finalized(&knights);
      if(!knights)
	MPI_Finalize();
    }
    int what_is_my_rank()
    {
      int rank;
      MPI_Comm_rank(FractalWorld,&rank);
      return rank;
    }
    int how_many_nodes()
    {
      int size;
      MPI_Comm_size(FractalWorld,&size);
      return size;
    }
    void FFTWStartup(const int& length_1,const bool& periodic)
    {
      cout << " FFTWStartup a " << FractalRank << endl;
      fftw_mpi_init();
      cout << " FFTWStartup b " << FractalRank << endl;
      bool wisdom_exists=fftw_import_wisdom_from_filename("wise/FractalWisdom.txt");
      if(periodic)
	{
	  const int length_c=(length_1+2)/2;
	  total_memory=fftw_mpi_local_size_3d(length_1,length_1,length_c,FractalWorld,&length_x,&start_x);
	  create_potRC();
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(length_1,length_1,length_1,potR,potC,FractalWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(length_1,length_1,length_1,potC,potR,FractalWorld,FFTW_MEASURE);
	}
      else
	{
	  const int length_11=length_1+1;
	  const int length_2=2*length_1;
	  const double g_c=pow(static_cast<double>(length_1),-5)/8.0;
	  cout << " g_c= " << g_c << " " << FractalRank << endl;
	  total_memory=fftw_mpi_local_size_3d(length_2,length_2,length_11,FractalWorld,&length_x,&start_x);
	  cout << " total_memory " << FractalRank << " " << total_memory << " " << length_x << " " << start_x << endl;
	  green.resize(length_x*length_11*length_11);
	  create_potRC();
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(length_2,length_2,length_2,potR,potC,FractalWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(length_2,length_2,length_2,potC,potR,FractalWorld,FFTW_MEASURE);
	  zeroR();
	  int length_22=length_2+2;
	  for(int nx=start_x;nx < start_x+length_x;++nx)
	    {
	      int nxa=min(nx,length_2-nx);
	      double x2=static_cast<double>(nxa*nxa)+0.25;
	      for(int ny=0;ny < length_2;++ny)
		{
		  int nya=min(ny,length_2-ny);
		  double y2=static_cast<double>(nya*nya);
		  for(int nz=0;nz<length_2;++nz)
		    {
		      int nza=min(nz,length_2-nz);
		      double z2=static_cast<double>(nza*nza);
		      double r2=z2+y2+x2;
		      potR[fftw_where(nx,ny,nz,length_2,length_22)]=-g_c/sqrt(r2);
		      //		      potR[fftw_where(nx,ny,nz,length_2,length_22)]=1.0;
		      //		      cout << " gra " << FractalRank << " " << nx <<  " " << ny <<  " " << nz << " " << g_c/sqrt(r2) << endl;
		    }	      
		}
	    }
	  //	  fftw_mpi_execute_dft_r2c(plan_rc,potR,potC);
	  fftw_execute(plan_rc);
	  for(int px=start_x;px<start_x+length_x;++px)
	    {
	      for(int py=0;py<length_11;++py)
		{
		  for(int pz=0;pz<length_11;++pz)
		    {
		      green[fftw_where(px,py,pz,length_11,length_11)]=potC[fftw_where(px,py,pz,length_2,length_11)][0];
		      //		      cout << " grb " << FractalRank << " " << px <<  " " << py <<  " " << pz << " ";
		      //		      cout << -green[fftw_where(px,py,pz,length_11,length_11)] << " ";
		      //		      cout << -potC[fftw_where(px,py,pz,length_2,length_11)][0] << " " << -potC[fftw_where(px,py,pz,length_2,length_11)][1] << endl;		    
		    }
		}
	    }
	  free_potRC();
	}
    }
    void FFTWFinal()
    {
      //      fftw_free(potR);
      //      fftw_free(potC);
      fftw_destroy_plan(plan_rc);
      fftw_destroy_plan(plan_cr);
      fftw_mpi_cleanup();
    }
    void dumpR(ofstream& FILE,const int& length,const bool& test)
    {
      int nx,ny,nz;
      for(nx=start_x;nx<start_x+length_x;nx++)
	{
	  for(ny=0;ny<length;ny++)
	    {
	      for(nz=0;nz<length;nz++)
		{
		  int n=fftw_where(nx,ny,nz,length,length+2);
		  if(!test || potR[n] != 0.0)
		    FILE << " dumpR " << nx << " " << ny << " " << nz << " " << n << " " << potR[n] << endl;
		}
	    }
	}
    }
    void create_potRC()
    {
      size_t sizeR=sizeof(double);
      size_t sizeC=sizeof(fftw_complex);
      potR=(double*) fftw_malloc(sizeR*2*total_memory);
      potC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
    }
    void free_potRC()
    {
      fftw_free(potR);
      fftw_free(potC);
    }
    void fftw_real_to_complex()
    {
      fftw_mpi_execute_dft_r2c(plan_rc,potR,potC);
    }
    void fftw_complex_to_real()
    {
      fftw_mpi_execute_dft_c2r(plan_cr,potC,potR);
    }
    inline int fftw_where(const int& i,const int& j,const int& k,const int& lb,const int& lc)
    {
      return k+(j+(i-start_x)*lb)*lc;
    }
    void calc_fftw_Slices(const int& length_a,const bool& periodic)
    {
      int paramsend[2]={start_x,start_x+length_x-1};
      int* paramrecv=(int*)malloc(2*FractalNodes*sizeof(int));
      int length_1=length_a;
      if(!periodic)
	length_1=2*length_a;
      paramsend[0]=start_x;
      paramsend[1]=start_x+length_x-1;
      cout << "calc_fftwa " << FractalRank << endl;
      MPI_Allgather(paramsend,2,MPI_INT,paramrecv,2,MPI_INT,FractalWorld);
      cout << "calc_fftwb " << FractalRank << endl;
      Slices.resize(FractalNodes);
      BoxS.resize(FractalNodes);
      BoxSL.resize(FractalNodes);
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  Slices[FR].resize(2);
	  Slices[FR][0]=paramrecv[2*FR];
	  Slices[FR][1]=paramrecv[2*FR+1];
	  BoxS[FR].resize(6);
	  BoxS[FR][0]=Slices[FR][0];
	  BoxS[FR][1]=Slices[FR][1];
	  BoxS[FR][2]=0;
	  BoxS[FR][3]=length_1-1;
	  BoxS[FR][4]=0;
	  BoxS[FR][5]=length_1-1;
	  BoxSL[FR].resize(3);
	  BoxSL[FR][0]=length_x;
	  BoxSL[FR][1]=length_1;
	  BoxSL[FR][2]=length_1;
	  cout << " slices " << FR << " " << Slices[FR][0] << " " << Slices[FR][1] << FractalRank << endl;
	}
      free(paramrecv);
      WhichSlice.resize(length_1);
      for(int nx=0;nx<length_1;nx++)
	{
	  bool success=false;
	  for(int S=0;S<FractalNodes;S++)
	    {
	      if(nx >= Slices[S][0] && nx <= Slices[S][1])
		{
		  WhichSlice[nx]=S;
		  success=true;
		  break;
		}
	    }
	  assert(success);
	}
      for(int ni=0;ni<length_1;ni++)
	cout << "whichslice " << FractalRank << " " << ni << " " << WhichSlice[ni] << endl;
    }
    void How_Many_On_Nodes(const int& count,vector <int>& counts)
    {
      counts.resize(FractalNodes);
      int* counts_out=new int[1];
      int* counts_in=new int[FractalNodes];
      counts_out[0]=count;
      MPI_Allgather(counts_out,1,MPI_INT,counts_in,1,MPI_INT,FractalWorld);
      for(int ni=0;ni<FractalNodes;ni++)
	counts[ni]=counts_in[ni];
      delete [] counts_out;
      delete [] counts_in;
    }
    void How_Many_Things_To_Send(vector <int>& counts_out_send,vector <int>& counts_in_send)
    {
      int* counts_out=new int[FractalNodes];
      int* counts_in=new int[FractalNodes];
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  counts_out[FR]=counts_out_send[FR];
	  cout << "sending a " << FractalRank << " " << FR << " " << counts_out[FR] << endl;
	}
      for(int FR=0;FR<FractalNodes;FR++)
	MPI_Gather(&counts_out[FR],1,MPI_INT,counts_in,1,MPI_INT,FR,FractalWorld);
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  counts_in_send[FR]=counts_in[FR];
	  cout << "sending b " << FractalRank << " " << FR << " " << counts_out[FR] << " " << counts_in[FR] << endl;
	}
      delete [] counts_out;
      delete [] counts_in;
    }
    void Send_Data_Somewhere(vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
			     vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
			     vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
    {
      int* displsI= new int[FractalNodes];
      int* displsR= new int[FractalNodes];
      int* countsI_in= new int[FractalNodes];
      int* countsI_out= new int[FractalNodes];
      int* countsR_in= new int[FractalNodes];
      int* countsR_out= new int[FractalNodes];
      displsI[0]=0;
      displsR[0]=0;
      int maxcount=-1;
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  countsI_in[FR]=counts_in_send[FR]*integers;
	  countsI_out[FR]=counts_out_send[FR]*integers;
	  countsR_in[FR]=counts_in_send[FR]*doubles;
	  countsR_out[FR]=counts_out_send[FR]*doubles;
	  maxcount=max(maxcount,counts_out_send[FR]);
	  if(FR > 0)
	    {
	      displsI[FR]=displsI[FR-1]+counts_in_send[FR-1]*integers;
	      displsR[FR]=displsR[FR-1]+counts_in_send[FR-1]*doubles;
	    }
	}
      how_manyI=displsI[FractalNodes-1]+counts_in_send[FractalNodes-1]*integers;
      how_manyR=displsR[FractalNodes-1]+counts_in_send[FractalNodes-1]*doubles;
      int* DataI_out=new int[max(maxcount*integers,1)];
      double* DataR_out=new double[max(maxcount*doubles,1)];
      int* dataI_in=new int[how_manyI];
      double* dataR_in=new double[how_manyR];
      cout << " howmanyIR " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  if(integers > 0)
	    {
	      for(int ni=0;ni<countsI_out[FR];ni++)
		DataI_out[ni]=dataI_out[FR][ni];
	      MPI_Gatherv(DataI_out,countsI_out[FR],MPI_INT,dataI_in,countsI_in,displsI,MPI_INT,FR,FractalWorld);  
	    }
	  if(doubles > 0)
	    {
	      for(int ni=0;ni<countsR_out[FR];ni++)
		DataR_out[ni]=dataR_out[FR][ni];
	      MPI_Gatherv(DataR_out,countsR_out[FR],MPI_DOUBLE,dataR_in,countsR_in,displsR,MPI_DOUBLE,FR,FractalWorld);  
	    }
	}
      delete [] displsI;
      delete [] displsR;
      delete [] countsI_in;
      delete [] countsI_out;
      delete [] countsR_in;
      delete [] countsR_out;
      delete [] DataI_out;
      delete [] DataR_out;
      displsI=0;
      displsR=0;
      countsI_in=0;
      countsI_out=0;
      countsR_in=0;
      countsR_out=0;
      DataI_out=0;
      DataR_out=0;
      cout << " how many " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
      if(integers > 0)
	{
	  dataI_in_send.resize(how_manyI);
	  for(int ni=0;ni<how_manyI;ni++)
	    dataI_in_send[ni]=dataI_in[ni];
	}
      if(doubles > 0)
	{
	  dataR_in_send.resize(how_manyR);
	  for(int ni=0;ni<how_manyR;ni++)
	    dataR_in_send[ni]=dataR_in[ni];
	}
      delete [] dataI_in;
      delete [] dataR_in;
      dataI_in=0;
      dataR_in=0;
    }
    void calc_total_particles(const int& NP)
    {
      int particles[1]={NP};
      Find_Sum_INT(particles,1);
      number_particles_total=particles[0];
    }
    void Find_Max_INT(int* integers,const int& how_long)
    {
      int* maxy=new int[how_long];
      MPI_Allreduce(integers,maxy,how_long,MPI_INT,MPI_MAX,FractalWorld);
      for(int ni=0;ni<how_long;ni++)
	integers[ni]=maxy[ni];
      delete [] maxy;
    }
    void Find_Max_DOUBLE(double* doubles,const int& how_long)
    {
      double* maxy=new double[how_long];
      MPI_Allreduce(doubles,maxy,how_long,MPI_DOUBLE,MPI_MAX,FractalWorld);
      for(int ni=0;ni<how_long;ni++)
	doubles[ni]=maxy[ni];
      delete [] maxy;
    }
    void Find_Sum_INT(int* integers,const int& how_long)
    {
      int* sumup=new int[how_long];
      MPI_Allreduce(integers,sumup,how_long,MPI_INT,MPI_SUM,FractalWorld);
      for(int ni=0;ni<how_long;ni++)
	integers[ni]=sumup[ni];
      delete [] sumup;
    }
    void Find_Sum_DOUBLE(vector <double>& doubles,const int& how_long)
    {
      double* doublesa=new double[how_long];
      for(int ni=0;ni<how_long;ni++)
	doublesa[ni]=doubles[ni];
      double* doublesb=new double[how_long];
      MPI_Allreduce(doublesa,doublesb,how_long,MPI_DOUBLE,MPI_SUM,FractalWorld);
      for(int ni=0;ni<how_long;ni++)
	doubles[ni]=doublesb[ni];
      delete [] doublesa;
      delete [] doublesb;
    }
    void Find_Sum_INT_to_Root(int* numbers,const int& how_long,const int& Root)
    {
      int* sumup=new int[how_long];
      MPI_Reduce(numbers,sumup,how_long,MPI_INT,MPI_SUM,Root,FractalWorld);
      for(int ni=0;ni<how_long;ni++)
	numbers[ni]=sumup[ni];
      delete [] sumup;
    }
    void Send_INT_from_Root(int* numbers,const int& how_long,const int& Root)
    {
      MPI_Bcast(numbers,how_long,MPI_INT,Root,FractalWorld);
    }
    void Full_Stop()
    {
      MPI_Barrier(HypreWorld);
    }
    void zeroR()
    {
      for(ptrdiff_t ni=0;ni<2*total_memory;ni++)
	potR[ni]=0.0;
    }
    double Clock()
    {
      return clock();
    }
  };
}
#endif
