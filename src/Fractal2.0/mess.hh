#ifndef _Mess_Defined_
#define _Mess_Defined_
namespace FractalSpace
{
  typedef ptrdiff_t pint;
  class Mess{
  public:
    int FractalRank;
    int FractalNodes;
    int FractalNodes0;
    int FractalNodes1;
    int FractalNodes2;
    int FFTRank;
    int FFTNodes;
    int HypreRank;
    int HypreNodes;
    int MPI_SWITCH;
    int MPI_MAX_COMMS;
    long int number_particles_total;
    Particle* parts_tmp;
    Particle* parts_tmpp;
    Particle* Parts_in;
    vector < vector <int> > Slices;
    vector < vector <int> > BoxS;
    vector < vector <int> > BoxSL;
    pint start_x;
    pint length_x;
    pint total_memory;
    fftw_plan plan_rc;
    fftw_plan plan_cr;
    double* potR;
    fftw_complex* potC; 
    vector <int> WhichSlice;
    /*
    vector <int>return_Slice_pos;
    vector <int>return_group;
    vector <int>return_point;
    vector <int>return_node;
    vector <int>what_Slice_point;
    */
    vector <double> green;
    MPI_Comm FractalWorld;
    MPI_Comm FFTWorld;
    MPI_Comm HypreWorld;
    MPI_Group FractalGroup;
    MPI_Group FFTGroup;
    MPI_Group HypreGroup;
    bool IAmAnFFTNode;
    bool IAmAHypreNode;
    vector <int>Hranks;
    vector <int>IHranks;
    vector <MPI_Comm> MComms;
    vector <MPI_Comm> HComms;
    vector <MPI_Group> HG;
    bool time_trial;
    bool standalone;
    File* p_file;
    double WallTime;
    double TreeTime;
    Mess():
      FractalRank(0),
      FractalNodes(1),
      FFTRank(0),
      FFTNodes(1234567),
      HypreRank(0),
      HypreNodes(0),
      MPI_SWITCH(512),
      MPI_MAX_COMMS(768),
      number_particles_total(-1),
      start_x(0),
      length_x(-1),
      total_memory(-1),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      standalone(true),
      TreeTime(-1.0)
    {
      WallTime=Clock();
      cout << " Empty Mess " << "\n";
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP,
	 int& FR0,int& FR1,int& FR2,const int& FN,MPI_Comm& FW):
      FractalRank(0),
      FractalNodes(1),
      FractalNodes0(FR0),
      FractalNodes1(FR1),
      FractalNodes2(FR2),
      FFTNodes(FN),
      HypreRank(0),
      HypreNodes(0),
      MPI_SWITCH(512),
      MPI_MAX_COMMS(768),
      number_particles_total(-1),
      start_x(0),
      length_x(GR),
      total_memory(1),
      FractalWorld(FW),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      TreeTime(-1.0)
    {
      //      cout << " Making a Mess with parameters" << "\n";
      int grid_length=GR;
      bool periodic=PR;
      WallTime=Clock();
      if(MR)
	{
	  MPIStartup(PR,FR0,FR1,FR2);
	  FractalRank=what_is_my_rank(); 
	  FractalNodes=how_many_nodes(); 
	  assert(FractalNodes == FR0*FR1*FR2);
	  FFTWStartup(grid_length,periodic);
	  calc_fftw_Slices(grid_length,periodic);	
	  calc_total_particles(NP);
	  make_MPI_Groups();
	}
      else
	{
	  number_particles_total=NP;
	  length_x=grid_length;
	}
      if(FractalRank == 0)
	cout << " made a mess " << FractalRank << " " << FractalNodes << " " << length_x << " " << start_x << " " << total_memory << "\n";
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP,const int& FN,MPI_Comm& FW):
      FractalRank(0),
      FractalNodes(1),
      FFTNodes(FN),
      HypreRank(0),
      HypreNodes(0),
      MPI_SWITCH(512),
      MPI_MAX_COMMS(768),
      number_particles_total(-1),
      start_x(0),
      length_x(GR),
      total_memory(1),
      FractalWorld(FW),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      TreeTime(-1.0)
    {
      //      cout << " Making a Mess with parameters" << "\n";
      int grid_length=GR;
      bool periodic=PR;
      WallTime=Clock();
      if(MR)
	{
	  MPIStartup();
	  FractalRank=what_is_my_rank(); 
	  FractalNodes=how_many_nodes(); 
	  FFTWStartup(grid_length,periodic);
	  calc_fftw_Slices(grid_length,periodic);	
	  calc_total_particles(NP);
	  make_MPI_Groups();
	}
      else
	{
	  number_particles_total=NP;
	  length_x=grid_length;
	}
      if(FractalRank == 0)
	cout << " made a mess " << FractalRank << " " << FractalNodes << " " << length_x << " " << start_x << " " << total_memory << "\n";
    }
    ~Mess()
    {
      cout << " starting to clean up a mess " << FractalRank << "\n";
      FFTWFinal();
      if(standalone)
	MPIFinal();
      cout << " cleaned up a mess " << FractalRank << "\n";
    }
    void MPIStartup()
    {
      //      cout << " Into MPIStartup " << "\n";
      int knights;
      MPI_Initialized(&knights);
      if(!knights)
	MPI_Init(NULL,NULL);
      FractalWorld=MPI_COMM_WORLD;
      FFTWorld=FractalWorld;
      HypreWorld=MPI_COMM_NULL;
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Comm_group(FFTWorld,&FFTGroup);
      //      MPI_Comm_group(HypreWorld,&HypreGroup);
      FractalRank=what_is_my_rank(); 
      FractalNodes=how_many_nodes();
      FFTRank=FractalRank;
      FFTNodes=min(FFTNodes,FractalNodes);
      HypreRank=FractalRank;
      HypreNodes=FractalNodes;
      if(FractalRank == 0)
	cout << " initialized MPI " << FractalRank << " " << FractalNodes << "\n";
    }
    void MPIStartup(bool PR,int& FR0,int& FR1,int& FR2)
    {
      //      int ranky;
      //      MPI_Comm_rank(FractalWorld,&ranky);
      //      cout << " Into MPIStartup A " << ranky << "\n";
      int knights;
      MPI_Initialized(&knights);
      if(!knights)
	MPI_Init(NULL,NULL);
      //      cout << " Into MPIStartup B " << ranky << "\n";
      int dims[]={FR0,FR1,FR2};
      int periods[]={PR,PR,PR};
      //      int periods[]={true,true,true};
      MPI_Cart_create(MPI_COMM_WORLD,3,dims,periods,true,&FractalWorld);
      FFTWorld=FractalWorld;
      HypreWorld=MPI_COMM_NULL;
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Comm_group(FFTWorld,&FFTGroup);
      FractalRank=what_is_my_rank(); 
      FractalNodes=how_many_nodes();
      FFTRank=FractalRank;
      FFTNodes=min(FFTNodes,FractalNodes);
      HypreRank=FractalRank;
      HypreNodes=FractalNodes;
      //      cout << " initialized MPI " << FractalRank << " " << FractalNodes << "\n";
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
    int what_is_my_rank(MPI_Comm& World)
    {
      int rank;
      MPI_Comm_rank(World,&rank);
      return rank;
    }
    int how_many_nodes(MPI_Comm& World)
    {
      int size;
      MPI_Comm_size(World,&size);
      return size;
    }
    int how_many_nodes()
    {
      int size;
      MPI_Comm_size(FractalWorld,&size);
      return size;
    }
    int what_is_my_FFT_rank()
    {
      int rank;
      if(FFTWorld == MPI_COMM_NULL)
	return -1;
      MPI_Comm_rank(FFTWorld,&rank);
      return rank;
    }
    int how_many_FFT_nodes()
    {
      int size;
      MPI_Comm_size(FFTWorld,&size);
      return size;
    }
    int what_is_my_Hypre_rank()
    {
      if(!IAmAHypreNode)
	return -1;
      int rank;
      MPI_Comm_rank(HypreWorld,&rank);
      return rank;
    }
    int how_many_Hypre_nodes()
    {
      if(!IAmAHypreNode)
	return -1;
      int size;
      MPI_Comm_size(HypreWorld,&size);
      return size;
    }
    void FFTWStartup(const int& length_1,const bool& periodic)
    {
      const pint Length_1=length_1;
      doFFTWorld(length_1,periodic);
      fftw_mpi_init();
      if(!IAmAnFFTNode) return;
      if(periodic)
	{
	  const pint Length_c=(Length_1+2)/2;
	  total_memory=fftw_mpi_local_size_3d(Length_1,Length_1,Length_c,FFTWorld,&length_x,&start_x);
	  cout << " total_memory " << FractalRank << " " << total_memory << " " << length_x << " " << start_x << "\n";
	  create_potRC();
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(Length_1,Length_1,Length_1,potR,potC,FFTWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(Length_1,Length_1,Length_1,potC,potR,FFTWorld,FFTW_MEASURE);
	  free_potRC();
	}
      else
	{
	  const pint Length_11=Length_1+1;
	  const pint Length_2=2*Length_1;
	  double g_c=pow(static_cast<double>(Length_1),-5)/8.0;
	  //	  cout << " g_c= " << g_c << " " << FractalRank << "\n";
	  total_memory=fftw_mpi_local_size_3d(Length_2,Length_2,Length_11,FFTWorld,&length_x,&start_x);
	  cout << " total_memory " << FractalRank << " " << total_memory << " " << length_x << " " << start_x << " " << g_c << "\n";
	  green.resize(length_x*Length_11*Length_11);
	  create_potRC();
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(Length_2,Length_2,Length_2,potR,potC,FFTWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(Length_2,Length_2,Length_2,potC,potR,FFTWorld,FFTW_MEASURE);
	  zeroR();
	  pint Length_22=Length_2+2;
	  for(pint nx=start_x;nx < start_x+length_x;++nx)
	    {
	      pint nxa=min(nx,Length_2-nx);
	      double x2=static_cast<double>(nxa*nxa)+0.25;
	      for(pint ny=0;ny < Length_2;++ny)
		{
		  pint nya=min(ny,Length_2-ny);
		  double y2=static_cast<double>(nya*nya);
		  for(pint nz=0;nz<Length_2;++nz)
		    {
		      pint nza=min(nz,Length_2-nz);
		      double z2=static_cast<double>(nza*nza);
		      double r2=z2+y2+x2;
		      potR[fftw_where(nx,ny,nz,Length_2,Length_22)]=-g_c/sqrt(r2);
		    }	      
		}
	    }
	  fftw_execute(plan_rc);
	  for(pint px=start_x;px<start_x+length_x;++px)
	    {
	      for(pint py=0;py<Length_11;++py)
		{
		  for(pint pz=0;pz<Length_11;++pz)
		    {
		      green[fftw_where(px,py,pz,Length_11,Length_11)]=potC[fftw_where(px,py,pz,Length_2,Length_11)][0];
		    }
		}
	    }
	  free_potRC();
	}
    }
    void doFFTWorld(int length_1,bool periodic)
    {
      //      cout << " messya " << FractalRank << " " << length_1 << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << IAmAnFFTNode << "\n";
      FFTRank=FractalRank;
      int howlong=length_1;
      if(!periodic)
	howlong*=2;
      int maxFFT=howlong/2;
      maxFFT=min(maxFFT,FractalNodes);
      FFTNodes=min(FFTNodes,maxFFT);
      while(howlong % FFTNodes != 0)
	{
	  FFTNodes--;
	  assert(FFTNodes > 1);
	}
      IAmAnFFTNode=FFTRank < FFTNodes;
      //      cout << " messyb " << FractalRank << " " << length_1 << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << maxFFT << " " << IAmAnFFTNode << "\n";
      vector <int>Franks(FFTNodes,-1);
      if(FFTNodes == FractalNodes)
	{
	  for(int ni=0;ni < FFTNodes;ni++)
	    Franks[ni]=ni;
	}
      else
	{
	  IAmAnFFTNode=false;
	  for(int ni=0;ni < FFTNodes;ni++)
	    {
	      Franks[ni]=ni;
	      IAmAnFFTNode=IAmAnFFTNode || Franks[ni] == FractalRank;
	    }
	}
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Group_incl(FractalGroup, FFTNodes,&(*Franks.begin()), &FFTGroup);
      MPI_Comm_create(FractalWorld, FFTGroup, &FFTWorld);
      if(!IAmAnFFTNode)
	{
	  FFTRank=-1;
	  start_x=9876543;
	  length_x=0;
	  total_memory=1;
	}
      else
	{
	  assert(FFTRank == what_is_my_FFT_rank());
	  assert(FFTNodes == how_many_FFT_nodes());
	}
      if(FractalRank == 0)
	cout << " messyc " << FractalRank << " " << length_1 << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << maxFFT << " " << IAmAnFFTNode << "\n";
    }
    void FFTWFinal()
    {
      if(IAmAnFFTNode)
	{
	  fftw_destroy_plan(plan_rc);
	  fftw_destroy_plan(plan_cr);
	}
      fftw_mpi_cleanup();
    }
    void dumpR(ofstream& FILE,const int& length)
    {
      pint nx,ny,nz;
      pint Length=length;
      for(nx=start_x;nx<start_x+length_x;nx++)
	{
	  for(ny=0;ny<Length;ny++)
	    {
	      for(nz=0;nz<Length;nz++)
		{
		  pint n=fftw_where(nx,ny,nz,length,length+2);
		  FILE << " dumpR " << nx << " " << ny << " " << nz << " " << n << " " << potR[n] << "\n";
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
    void create_potR()
    {
      size_t sizeR=sizeof(double);
      potR=(double*) fftw_malloc(sizeR*2*total_memory);
    }
    void free_potR()
    {
      fftw_free(potR);
    }
    void create_potC()
    {
      size_t sizeC=sizeof(fftw_complex);
      potC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
    }
    void free_potC()
    {
      fftw_free(potC);
    }
    void fftw_real_to_complex()
    {
      if(IAmAnFFTNode)
	fftw_mpi_execute_dft_r2c(plan_rc,potR,potC);
    }
    void fftw_complex_to_real()
    {
      if(IAmAnFFTNode)
	fftw_mpi_execute_dft_c2r(plan_cr,potC,potR);
    }
    inline int fftw_where(const int& i,const int& j,const int& k,const int& lb,const int& lc)
    {
      return k+(j+(i-start_x)*lb)*lc;
    }
    void calc_fftw_Slices(const int& length_a,const bool& periodic)
    {
      //      int paramsend[2]={(int)start_x,(int)(start_x+length_x-1)};
      //      int* paramrecv=new int[2*FractalNodes];
      vector <int>paramsend(2);
      paramsend[0]=start_x;
      paramsend[1]=start_x+length_x-1;
      vector <int>paramrecv(2*FractalNodes);
      int length_1=length_a;
      if(!periodic)
	length_1=2*length_a;
      //      cout << "calc_fftwa " << FFTRank << " " << start_x << " " << length_x << "\n";
      my_AllgatherI(paramsend,paramrecv,2);
      //      MPI_Allgather(paramsend,2,MPI_INT,paramrecv,2,MPI_INT,FractalWorld);
      cout << "calc_fftwb " << FFTRank << " " << start_x << " " << length_x << "\n";
      Slices.resize(FractalNodes); // this is not an error.
      BoxS.resize(FractalNodes); // it must be dimensioned
      BoxSL.resize(FractalNodes); // this way
      for(int FR=0;FR<FFTNodes;FR++)
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
	  if(FFTRank == 0)
	    cout << " slices " << FFTRank << " " << Slices[FR][0] << " " << Slices[FR][1] << " " << FR << " " << FractalRank << "\n";
	}
      WhichSlice.assign(length_1,-10);
      bool allok=true;
      for(int nx=0;nx<length_1;nx++)
	{
	  bool success=false;
	  for(int S=0;S<FFTNodes;S++)
	    {
	      if(nx >= Slices[S][0] && nx <= Slices[S][1])
		{
		  WhichSlice[nx]=S;
		  success=true;
		  break;
		}
	    }
	  if(!success)
	    {
	      allok=false;
	      for(int nx=0;nx<length_1;nx++)
		if(FFTRank == 0) cout << " success " << FractalRank << " " << FFTRank << " " << nx << " " << WhichSlice[nx] << "\n";
	    }
	}
      for(int ni=0;ni<length_1;ni++)
	if(FFTRank == 0) cout << "whichslice " << FFTRank << " " << ni << " " << WhichSlice[ni] << "\n";
      assert(allok);
    }
    template <class T> void How_Many_On_Nodes(T count,vector <T>& counts)
    {
      counts.resize(FractalNodes);
      vector <T>counts_out;
      counts_out.push_back(count);
      my_AllgatherI(counts_out,counts,1);
    }
    void How_Many_Things_To_Send_I(vector <int>& counts_out_send,vector <int>& counts_in_send)
    {
      How_Many_Things_To_Send_I(FractalWorld,counts_out_send,counts_in_send);
    }
    void How_Many_Things_To_Send_I(MPI_Comm& World,
				   vector <int>& counts_out_send,vector <int>& counts_in_send)
    {
      int Nodes;
      MPI_Comm_size(World,&Nodes);
      vector <vector <int> > dataI_out(FractalNodes);
      vector <vector <double> > dataR_out(FractalNodes);
      vector <double> dataR_in;
      vector <int> dataI_in;
      int how_manyI=-1;
      int how_manyR=-1;
      vector <int>counts_out(Nodes,1);
      vector <int> counts_in(Nodes,1);
      for(int FR=0;FR<Nodes;FR++)
	dataI_out[FR].push_back(counts_out_send[FR]);
      Send_Data_Somewhere_No_Block(World,counts_out,counts_in,1,0,dataI_out,dataI_in,how_manyI,dataR_out,dataR_in,how_manyR);
      counts_in_send=dataI_in;
    }
    void Send_Data_Somewhere_No_Block(vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
				      vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
				      vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
    {
      ofstream& FF=p_file->FileFractal;
      vector <int>displsI(FractalNodes,0);
      vector <int>displsR(FractalNodes,0);
      vector <int>countsI_in(FractalNodes,0);
      vector <int>countsR_in(FractalNodes,0);
      vector <int>countsI_out(FractalNodes,0);
      vector <int>countsR_out(FractalNodes,0);
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  countsI_in[FR]=counts_in_send[FR]*integers;
	  countsI_out[FR]=counts_out_send[FR]*integers;
	  countsR_in[FR]=counts_in_send[FR]*doubles;
	  countsR_out[FR]=counts_out_send[FR]*doubles;
	  if(FR > 0)
	    {
	      displsI[FR]=displsI[FR-1]+counts_in_send[FR-1]*integers;
	      displsR[FR]=displsR[FR-1]+counts_in_send[FR-1]*doubles;
	    }
	}
      how_manyI=displsI[FractalNodes-1]+counts_in_send[FractalNodes-1]*integers;
      how_manyR=displsR[FractalNodes-1]+counts_in_send[FractalNodes-1]*doubles;
      int extraI=countsI_out[FractalRank];
      int extraR=countsR_out[FractalRank];
      int startI=displsI[FractalRank];
      int startR=displsR[FractalRank];
      const int tagI=0;
      const int tagR=1;
      vector <MPI_Request> requestIout;
      vector <MPI_Request> requestRout;
      vector <MPI_Request> requestIin;
      vector <MPI_Request> requestRin;
      FF << " howmanyIR " << FractalRank << " " << how_manyI << " " << how_manyR << "\n";
      int answer=-1;
      if(integers > 0)
	{
	  dataI_in_send.resize(how_manyI);
	  for(int ni=0;ni<extraI;ni++)
	    dataI_in_send[ni+startI]=dataI_out[FractalRank][ni];
	  //	  std::copy(dataI_out[FractalRank].begin(),dataI_out[FractalRank].begin()+extraI,dataI_in_send.begin()+startI)
	}
      if(doubles > 0)
	{
	  dataR_in_send.resize(how_manyR);
	  for(int ni=0;ni<extraR;ni++)
	    dataR_in_send[ni+startR]=dataR_out[FractalRank][ni];
	  //	  std::copy(dataR_out[FractalRank].begin(),dataR_out[FractalRank].begin()+extraR,dataR_in_send.begin()+startR)
	}
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  if(FR == FractalRank || counts_in_send[FR] == 0)
	    continue;
	  if(integers > 0)
	    {
	      requestIin.push_back(MPI_Request());
	      answer=MPI_Irecv(&(*(dataI_in_send.begin()+displsI[FR])),countsI_in[FR],MPI_INT,FR,tagI,FractalWorld,&requestIin.back());
	      MPI_MYTest(0,answer);
	    }
	  if(doubles > 0)
	    {
	      requestRin.push_back(MPI_Request());
	      answer=MPI_Irecv(&(*(dataR_in_send.begin()+displsR[FR])),countsR_in[FR],MPI_DOUBLE,FR,tagR,FractalWorld,&requestRin.back());
	      MPI_MYTest(1,answer);
	    }
	}
      //
      Full_Stop_Do_Not_Argue();
      //
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  if(FR == FractalRank || counts_out_send[FR] == 0)
	    continue;
	  if(integers > 0)
	    {
	      requestIout.push_back(MPI_Request());
	      answer=MPI_Isend(&(*dataI_out[FR].begin()),countsI_out[FR],MPI_INT,FR,tagI,FractalWorld,&requestIout.back());
	      MPI_MYTest(2,answer);
	    }
	  if(doubles > 0)
	    {
	      requestRout.push_back(MPI_Request());
	      answer=MPI_Isend(&(*dataR_out[FR].begin()),countsR_out[FR],MPI_DOUBLE,FR,tagR,FractalWorld,&requestRout.back());
	      MPI_MYTest(3,answer);
	    }
	}
      if(integers > 0)
	{
	  vector <MPI_Status> statusIout(requestIout.size());
	  vector <MPI_Status> statusIin(requestIin.size());
	  answer=MPI_Waitall(requestIout.size(),&(*requestIout.begin()),&(*statusIout.begin()));
	  MPI_MYTest(4,answer);
	  answer=MPI_Waitall(requestIin.size(),&(*requestIin.begin()),&(*statusIin.begin()));
	  MPI_MYTest(5,answer);
	}
      if(doubles > 0)
	{
	  vector <MPI_Status> statusRout(requestRout.size());
	  vector <MPI_Status> statusRin(requestRin.size());
	  answer=MPI_Waitall(requestRout.size(),&(*requestRout.begin()),&(*statusRout.begin()));
	  MPI_MYTest(6,answer);
	  answer=MPI_Waitall(requestRin.size(),&(*requestRin.begin()),&(*statusRin.begin()));
	  MPI_MYTest(7,answer);
	}
      FF << " how many " << FractalRank << " " << how_manyI << " " << how_manyR << "\n";
    }
    void Send_Data_Somewhere_No_Block(MPI_Comm& World,
				      vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
				      vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
				      vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
    {
      int Rank;
      MPI_Comm_rank(World,&Rank);
      int Nodes;
      MPI_Comm_size(World,&Nodes);
      ofstream& FF=p_file->FileFractal;
      vector <int>displsI(Nodes,0);
      vector <int>displsR(Nodes,0);
      vector <int>countsI_in(Nodes,0);
      vector <int>countsR_in(Nodes,0);
      vector <int>countsI_out(Nodes,0);
      vector <int>countsR_out(Nodes,0);
      for(int FR=0;FR<Nodes;FR++)
	{
	  countsI_in[FR]=counts_in_send[FR]*integers;
	  countsI_out[FR]=counts_out_send[FR]*integers;
	  countsR_in[FR]=counts_in_send[FR]*doubles;
	  countsR_out[FR]=counts_out_send[FR]*doubles;
	  if(FR > 0)
	    {
	      displsI[FR]=displsI[FR-1]+counts_in_send[FR-1]*integers;
	      displsR[FR]=displsR[FR-1]+counts_in_send[FR-1]*doubles;
	    }
	}
      how_manyI=displsI[Nodes-1]+counts_in_send[Nodes-1]*integers;
      how_manyR=displsR[Nodes-1]+counts_in_send[Nodes-1]*doubles;
      int extraI=countsI_out[Rank];
      int extraR=countsR_out[Rank];
      int startI=displsI[Rank];
      int startR=displsR[Rank];
      const int tagI=0;
      const int tagR=1;
      vector <MPI_Request> requestIout;
      vector <MPI_Request> requestRout;
      vector <MPI_Request> requestIin;
      vector <MPI_Request> requestRin;
      FF << " howmanyIR " << Rank << " " << how_manyI << " " << how_manyR << "\n";
      int answer=-1;
      if(integers > 0)
	{
	  dataI_in_send.resize(how_manyI);
	  for(int ni=0;ni<extraI;ni++)
	    dataI_in_send[ni+startI]=dataI_out[Rank][ni];
	  //	  std::copy(dataI_out[Rank].begin(),dataI_out[Rank].begin()+extraI,dataI_in_send.begin()+startI)
	}
      if(doubles > 0)
	{
	  dataR_in_send.resize(how_manyR);
	  for(int ni=0;ni<extraR;ni++)
	    dataR_in_send[ni+startR]=dataR_out[Rank][ni];
	  //	  std::copy(dataR_out[Rank].begin(),dataR_out[Rank].begin()+extraR,dataR_in_send.begin()+startR)
	}
      for(int FR=0;FR<Nodes;FR++)
	{
	  if(FR == Rank || counts_in_send[FR] == 0)
	    continue;
	  if(integers > 0)
	    {
	      requestIin.push_back(MPI_Request());
	      answer=MPI_Irecv(&(*(dataI_in_send.begin()+displsI[FR])),countsI_in[FR],MPI_INT,FR,tagI,World,&requestIin.back());
	      MPI_MYTest(0,answer);
	    }
	  if(doubles > 0)
	    {
	      requestRin.push_back(MPI_Request());
	      answer=MPI_Irecv(&(*(dataR_in_send.begin()+displsR[FR])),countsR_in[FR],MPI_DOUBLE,FR,tagR,World,&requestRin.back());
	      MPI_MYTest(1,answer);
	    }
	}
      //
      Full_Stop_Do_Not_Argue(World);
      //
      for(int FR=0;FR<Nodes;FR++)
	{
	  if(FR == Rank || counts_out_send[FR] == 0)
	    continue;
	  if(integers > 0)
	    {
	      requestIout.push_back(MPI_Request());
	      answer=MPI_Isend(&(*dataI_out[FR].begin()),countsI_out[FR],MPI_INT,FR,tagI,World,&requestIout.back());
	      MPI_MYTest(2,answer);
	    }
	  if(doubles > 0)
	    {
	      requestRout.push_back(MPI_Request());
	      answer=MPI_Isend(&(*dataR_out[FR].begin()),countsR_out[FR],MPI_DOUBLE,FR,tagR,World,&requestRout.back());
	      MPI_MYTest(3,answer);
	    }
	}
      if(integers > 0)
	{
	  vector <MPI_Status> statusIout(requestIout.size());
	  vector <MPI_Status> statusIin(requestIin.size());
	  answer=MPI_Waitall(requestIout.size(),&(*requestIout.begin()),&(*statusIout.begin()));
	  MPI_MYTest(4,answer);
	  answer=MPI_Waitall(requestIin.size(),&(*requestIin.begin()),&(*statusIin.begin()));
	  MPI_MYTest(5,answer);
	}
      if(doubles > 0)
	{
	  vector <MPI_Status> statusRout(requestRout.size());
	  vector <MPI_Status> statusRin(requestRin.size());
	  answer=MPI_Waitall(requestRout.size(),&(*requestRout.begin()),&(*statusRout.begin()));
	  MPI_MYTest(6,answer);
	  answer=MPI_Waitall(requestRin.size(),&(*requestRin.begin()),&(*statusRin.begin()));
	  MPI_MYTest(7,answer);
	}
      FF << " how many " << Rank << " " << how_manyI << " " << how_manyR << "\n";
    }
    void make_MPI_Hypre_Groups()
    {
      //      cout << " making hypre " << FractalRank << " " << HypreRank << "\n";
      HG.resize(3);
      HComms.clear();
      HComms.resize(3);
      double aNodes=HypreNodes;
      int HypreNodes0=pow(aNodes-0.5,1.0/3.0)+1.0;
      double a12=HypreNodes/HypreNodes0;
      int HypreNodes1=sqrt(a12-0.5)+1.0;
      int HypreNodes01=HypreNodes0*HypreNodes1;
      int HypreNodes2=HypreNodes/HypreNodes01;
      int HypreRank0=HypreRank % HypreNodes0;
      int HypreRank1=(HypreRank/HypreNodes0) % HypreNodes1;
      int HypreRank2=HypreRank/HypreNodes01;
      int HypreNodesBox=HypreNodes01*HypreNodes2;
      int extras=HypreNodes-HypreNodesBox;
      int ExtraLines=extras/HypreNodes0;
      int ExtraNodes=extras % HypreNodes0;
      if(HypreRank == 0)
	{
	  cout << " MAKE HGD " << HypreNodes << " " << HypreNodes0 << " " << HypreNodes1 << " " << HypreNodes2 << " " << HypreNodes01;
	  cout << " " << HypreRank0 <<  " " << HypreRank1 <<  " " << HypreRank2 << " " << ExtraLines << " " << ExtraNodes << "\n";
	}
      vector < vector <int> > RanksH;
      RanksH.clear();
      RanksH.resize(3);
      int rr=HypreRank % HypreNodes01;
      while(rr < HypreNodes)
	{
	  RanksH[2].push_back(rr);
	  rr+=HypreNodes01;
	}
      if(HypreRank2 < HypreNodes2)
	{
	  int ni=0;
	  int rr=HypreRank0+HypreRank2*HypreNodes01;
	  while(ni < HypreNodes1)
	    {
	      RanksH[1].push_back(rr);
	      rr+=HypreNodes0;
	      ni++;
	    }
	}
      if(HypreRank2 < HypreNodes2 || HypreRank1 < ExtraLines)
	{
	  int ni=0;
	  rr=HypreRank1*HypreNodes0+HypreRank2*HypreNodes01;
	  while(ni < HypreNodes0)
	    {
	      RanksH[0].push_back(rr);
	      rr++;
	      ni++;
	    }
	}
      if(HypreRank2 == HypreNodes2)
	{
	  int rr=HypreRank0+HypreRank2*HypreNodes01;
	  int ni=0;
	  while(ni < ExtraLines)
	    {
	      RanksH[1].push_back(rr);
	      rr+=HypreNodes0;
	      ni++;
	    }
	  if(HypreRank0 < ExtraNodes)
	    RanksH[1].push_back(rr);
	  if(ExtraNodes > 0 && ((HypreRank1 == ExtraLines && ExtraLines > 0)|| ExtraLines == 0)) 
	    {
	      int ni=0;
	      rr=HypreRank1*HypreNodes0+HypreRank2*HypreNodes01;
	      while(ni < ExtraNodes)
		{
		  RanksH[0].push_back(rr);
		  rr++;
		  ni++;
		}
	    }
	}
      MPI_Group_incl(HypreGroup,RanksH[2].size(),&(*RanksH[2].begin()),&HG[2]);
      MPI_Comm_create(HypreWorld,HG[2],&HComms[2]);
      MPI_Group_incl(HypreGroup,RanksH[1].size(),&(*RanksH[1].begin()),&HG[1]);
      MPI_Comm_create(HypreWorld,HG[1],&HComms[1]);
      MPI_Group_incl(HypreGroup,RanksH[0].size(),&(*RanksH[0].begin()),&HG[0]);
      MPI_Comm_create(HypreWorld,HG[0],&HComms[0]);
    }
    void make_MPI_Groups()
    {
      if(FractalNodes <= MPI_SWITCH)
	return;
      vector <MPI_Group> MG;
      int FractalNodes01=FractalNodes0*FractalNodes1;
      int FractalRank0=FractalRank % FractalNodes0;
      int FractalRank1=(FractalRank/FractalNodes0) % FractalNodes1;
      int FractalRank2=FractalRank/FractalNodes01;
      vector <int>Ranks0(FractalNodes0);
      vector <int>Ranks1(FractalNodes1);
      vector <int>Ranks2(FractalNodes2);

      for(int FR0=0;FR0<FractalNodes0;FR0++)
	Ranks0[FR0]=FR0+(FractalRank1+FractalRank2*FractalNodes1)*FractalNodes0;
      MComms.push_back(MPI_Comm());
      MG.push_back(MPI_Group());
      MPI_Group_incl(FractalGroup,FractalNodes0,&(*Ranks0.begin()),&MG.back());
      MPI_Comm_create(FractalWorld,MG.back(),&MComms.back());

      for(int FR1=0;FR1<FractalNodes1;FR1++)
	Ranks1[FR1]=FractalRank0+(FR1+FractalRank2*FractalNodes1)*FractalNodes0;
      MComms.push_back(MPI_Comm());
      MG.push_back(MPI_Group());
      MPI_Group_incl(FractalGroup,FractalNodes1,&(*Ranks1.begin()),&MG.back());
      MPI_Comm_create(FractalWorld,MG.back(),&MComms.back());

      for(int FR2=0;FR2<FractalNodes2;FR2++)
	Ranks2[FR2]=FractalRank0+(FractalRank1+FR2*FractalNodes1)*FractalNodes0;

      MComms.push_back(MPI_Comm());
      MG.push_back(MPI_Group());
      MPI_Group_incl(FractalGroup,FractalNodes2,&(*Ranks2.begin()),&MG.back());
      MPI_Comm_create(FractalWorld,MG.back(),&MComms.back());
    }
    void Send_Data_Some_How(int tag,vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
			    vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
			    vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
    {
      Send_Data_Some_How(tag,FractalWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in_send,how_manyI,
				   dataR_out,dataR_in_send,how_manyR);
    }
    void Send_Data_Some_How(int tag,MPI_Comm& World,
			    vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
			    vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
			    vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
    {
      int Rank=what_is_my_rank(World);
      int Nodes=how_many_nodes(World);
      bool small=Nodes <= MPI_SWITCH;
      bool foreign=World != FractalWorld && World != HypreWorld;
      p_file->DUMPS << " tag= " << tag << " " << Rank << " " << Nodes << " " << small << " " << foreign << endl;
      cout << " tag= " << tag << " " << Rank << " " << Nodes << " " << small << " " << foreign << endl;
      if(Rank == 0)
	cout << " SOMEWHOW " << FractalRank << " " << Nodes << " " << MPI_SWITCH << " " << small << foreign << " ";
      if(small || foreign)
	{
	  if(Rank == 0)
	    cout << "A" << "\n";
	  How_Many_Things_To_Send_I(World,counts_out,counts_in);
	  Send_Data_Somewhere_No_Block(World,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in_send,how_manyI,
				       dataR_out,dataR_in_send,how_manyR);
	  if(Rank == 0)
	    cout << "AA" << "\n";
	}
      else if(World == HypreWorld)
	{
	  if(Rank == 0)
	    cout << "B" << "\n";
	  Send_Data_Hypre_Directions(counts_out,counts_in,integers,doubles,
				     dataI_out,dataI_in_send,how_manyI,
				     dataR_out,dataR_in_send,how_manyR);
	  if(Rank == 0)
	    cout << "BB" << "\n";
	}
      else if(tag==0 || tag==4)
	{
	  if(Rank == 0)
	    cout << "C" << "\n";
	  Send_Data_Other_Directions(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in_send,how_manyI,
				   dataR_out,dataR_in_send,how_manyR);
	  if(Rank == 0)
	    cout << "CC" << "\n";
	}
      else
	{
	  if(Rank == 0)
	    cout << "D" << "\n";
	  Send_Data_One_Directions(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in_send,how_manyI,
				   dataR_out,dataR_in_send,how_manyR);
	  if(Rank == 0)
	    cout << "DD" << "\n";
	}
    }
    //
    void Send_Data_One_Directions(vector <int>& counts_out,vector <int>& counts_in,int integers,int doubles,
				  vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
				  vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR)
    {
      ofstream& FF=p_file->DUMPS;
      int FractalNodes01=FractalNodes0*FractalNodes1;
      int FractalRank2=FractalRank/FractalNodes01;
      dataI_in.clear();
      dataR_in.clear();
      vector <int>countsa_out(FractalNodes2,0);
      vector <int>countsa_in(FractalNodes2);
      FF << " Send AA " << FractalRank << "\n";
      int totals=0;
      int nIdata=0;
      int nRdata=0;
      try
	{
	  vector <int> tmpI0=dataI_out[0];
	  dataI_out[0].clear();
	  nIdata=0;
	  countsa_out[0]+=counts_out[0];
	  for(int ni=0;ni<counts_out[0];ni++)
	    {
	      dataI_out[0].push_back(0);
	      for(int ints=0;ints<integers;ints++)
		{
		  dataI_out[0].push_back(tmpI0[nIdata]);
		  nIdata++;
		}
	    }
	  tmpI0.clear();
	  for(int FR=1;FR<FractalNodes;FR++)
	    {
	      int FR2=FR/FractalNodes01;
	      countsa_out[FR2]+=counts_out[FR];
	      nIdata=0;
	      nRdata=0;
	      for(int ni=0;ni<counts_out[FR];ni++)
		{
		  dataI_out[FR2].push_back(FR);
		  for(int ints=0;ints<integers;ints++)
		    {
		      dataI_out[FR2].push_back(dataI_out[FR][nIdata]);
		      nIdata++;
		    }
		  for(int reals=0;reals<doubles;reals++)
		    {
		      dataR_out[FR2].push_back(dataR_out[FR][nRdata]);
		      nRdata++;
		    }
		  totals++;
		}
	      dataI_out[FR].clear();
	      dataR_out[FR].clear();
	    }
	  dataI_out.resize(FractalNodes2);
	  dataR_out.resize(FractalNodes2);
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT A " << FractalRank << " " << ba.what() << " " << totals << "\n";
	  FF << " DUMP IT A " << ba.what() << " " << totals << "\n";
	  for(int FR2=0;FR2<FractalNodes2;FR2++)
	    FF << FR2 << " " << dataI_out[FR2].size() << " " << dataR_out[FR2].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      FF << " Send BB " << FractalRank << "\n";
      How_Many_Things_To_Send_I(MComms[2],countsa_out,countsa_in);
      int total_in=0;
      int total_out=0;
      for(int FR2=0;FR2<FractalNodes2;FR2++)
	{
	  total_out+=countsa_out[FR2];
	  total_in+=countsa_in[FR2];
	}
      FF << " TotalsOne 2 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << endl;
      Send_Data_Somewhere_No_Block(MComms[2],countsa_out,countsa_in,
				   integers+1,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
      FF << " Send CC " << FractalRank << endl;
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes1);
      dataR_out.resize(FractalNodes1);
      countsa_out.assign(FractalNodes1,0);
      int countI=0;
      int countR=0;
      try
	{
	  for(int FR2=0;FR2<FractalNodes2;FR2++)
	    {
	      int FRFrom=FractalRank+(FR2-FractalRank2)*FractalNodes01;
	      for(int c=0;c<countsa_in[FR2];c++)
		{
		  int FR=dataI_in[countI];
		  countI++;
		  int FR1=(FR/FractalNodes0) % FractalNodes1;
		  countsa_out[FR1]++;
		  dataI_out[FR1].push_back(FR);
		  dataI_out[FR1].push_back(FRFrom);
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_out[FR1].push_back(dataI_in[countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_out[FR1].push_back(dataR_in[countR]);
		      countR++;
		    }
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT C " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT C " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes2;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      Full_Stop_Do_Not_Argue(MComms[2]);
      Full_Stop_Do_Not_Argue(MComms[1]);
      FF << " Send DD " << FractalRank << "\n";
      countsa_in.assign(FractalNodes1,0);
      How_Many_Things_To_Send_I(MComms[1],countsa_out,countsa_in);
      dataI_in.clear();
      dataR_in.clear();

      total_in=0;
      total_out=0;
      for(int FR1=0;FR1<FractalNodes1;FR1++)
	{
	  total_out+=countsa_out[FR1];
	  total_in+=countsa_in[FR1];
	}
      FF << " TotalsOne 1 " << total_out*(integers+2) << " " << total_out*doubles << " " <<  total_in*(integers+2) << " " << total_in*doubles << endl;
      Send_Data_Somewhere_No_Block(MComms[1],countsa_out,countsa_in,
				   integers+2,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);

      FF << " Send EE " << FractalRank << endl;
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes0);
      dataR_out.resize(FractalNodes0);
      countsa_out.assign(FractalNodes0,0);
      countI=0;
      countR=0;
      try
	{
	  for(int FR1=0;FR1<FractalNodes1;FR1++)
	    {
	      for(int c=0;c<countsa_in[FR1];c++)
		{
		  int FR=dataI_in[countI];
		  countI++;
		  int FR0=FR % FractalNodes0;
		  countsa_out[FR0]++;
		  int FRFrom=dataI_in[countI];
		  countI++;
		  dataI_out[FR0].push_back(FRFrom);
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_out[FR0].push_back(dataI_in[countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_out[FR0].push_back(dataR_in[countR]);
		      countR++;
		    }
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT D " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT D " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes1;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      Full_Stop_Do_Not_Argue(MComms[1]);
      Full_Stop_Do_Not_Argue(MComms[0]);
      FF << " Send FF " << FractalRank << "\n";
      countsa_in.assign(FractalNodes0,0);
      How_Many_Things_To_Send_I(MComms[0],countsa_out,countsa_in);
      dataI_in.clear();
      dataR_in.clear();
      total_in=0;
      total_out=0;
      for(int FR0=0;FR0<FractalNodes0;FR0++)
	{
	  total_out+=countsa_out[FR0];
	  total_in+=countsa_in[FR0];
	}
      FF << " TotalsOne 0 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << endl;
      Send_Data_Somewhere_No_Block(MComms[0],countsa_out,countsa_in,
				   integers+1,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);

      FF << " Send GG " << FractalRank << endl;
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes);
      dataR_out.resize(FractalNodes);
      counts_in.assign(FractalNodes,0);
      countI=0;
      countR=0;
      try
	{
	  for(int FR0=0;FR0<FractalNodes0;FR0++)
	    {
	      for(int c=0;c<countsa_in[FR0];c++)
		{
		  int FRFrom=dataI_in[countI];
		  countI++;
		  counts_in[FRFrom]++;
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_out[FRFrom].push_back(dataI_in[countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_out[FRFrom].push_back(dataR_in[countR]);
		      countR++;
		    }
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT E " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT E " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes0;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      dataI_in.clear();
      dataR_in.clear();
      FF << " Send HH " << FractalRank <<  "\n";
      how_manyI=0;
      how_manyR=0;
      try
	{
	  for(int FR=0;FR<FractalNodes;FR++)
	    {
	      countI=0;
	      countR=0;
	      for(int c=0;c<counts_in[FR];c++)
		{
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_in.push_back(dataI_out[FR][countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_in.push_back(dataR_out[FR][countR]);
		      countR++;
		    }
		}
	      how_manyI+=countI;
	      how_manyR+=countR;
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT F " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT F " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      Full_Stop_Do_Not_Argue();
      FF << " Send II " << FractalRank << " " << how_manyI << " " << how_manyR << "\n";
    }
    //
    void Send_Data_Other_Directions(vector <int>& counts_out,vector <int>& counts_in,int integers,int doubles,
				  vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
				  vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR)
    {
      ofstream& FF=p_file->DUMPS;
      int FractalNodes01=FractalNodes0*FractalNodes1;
      int FractalRank0=FractalRank % FractalNodes0;
      dataI_in.clear();
      dataR_in.clear();
      vector <int>countsa_out(FractalNodes0,0);
      vector <int>countsa_in(FractalNodes0);
      FF << " Send AA " << FractalRank << "\n";
      int totals=0;
      try
	{
	  for(int FR=0;FR<FractalNodes;FR++)
	    {
	      int FR0=FR % FractalNodes0;;
	      countsa_out[FR0]+=counts_out[FR];
	      int nIdata=0;
	      int nRdata=0;
	      for(int ni=0;ni<counts_out[FR];ni++)
		{
		  dataI_in.push_back(FR);
		  for(int ints=0;ints<integers;ints++)
		    {
		      dataI_in.push_back(dataI_out[FR][nIdata]);
		      nIdata++;
		    }
		  for(int reals=0;reals<doubles;reals++)
		    {
		      dataR_in.push_back(dataR_out[FR][nRdata]);
		      nRdata++;
		    }
		  totals++;
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT A " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT A " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes0);
      dataR_out.resize(FractalNodes0);
      int counterI=0;
      int counterR=0;
      try
	{
	  for(int ni=0;ni<totals;ni++)
	    {
	      int FR=dataI_in[counterI];
	      int FR0=FR % FractalNodes0;
	      counterI++;
	      dataI_out[FR0].push_back(FR);
	      for(int niI=0;niI<integers;niI++)
		{
		  dataI_out[FR0].push_back(dataI_in[counterI]);
		  counterI++;
		}
	      for(int niR=0;niR<doubles;niR++)
		{
		  dataR_out[FR0].push_back(dataR_in[counterR]);
		  counterR++;
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT B " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT B " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes0;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      dataI_in.clear();
      dataR_in.clear();
      FF << " Send BB " << FractalRank << "\n";
      How_Many_Things_To_Send_I(MComms[0],countsa_out,countsa_in);

      int total_in=0;
      int total_out=0;
      for(int FR0=0;FR0<FractalNodes0;FR0++)
	{
	  total_out+=countsa_out[FR0];
	  total_in+=countsa_in[FR0];
	}
      FF << " TotalsOther 0 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << endl;

      Send_Data_Somewhere_No_Block(MComms[0],countsa_out,countsa_in,
				   integers+1,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
      FF << " Send CC " << FractalRank << endl;
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes1);
      dataR_out.resize(FractalNodes1);
      countsa_out.assign(FractalNodes1,0);
      int countI=0;
      int countR=0;
      try
	{
	  for(int FR0=0;FR0<FractalNodes0;FR0++)
	    {
	      int FRFrom=FractalRank+(FR0-FractalRank0);
	      for(int c=0;c<countsa_in[FR0];c++)
		{
		  int FR=dataI_in[countI];
		  countI++;
		  int FR1=(FR/FractalNodes0) % FractalNodes1;
		  countsa_out[FR1]++;
		  dataI_out[FR1].push_back(FR);
		  dataI_out[FR1].push_back(FRFrom);
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_out[FR1].push_back(dataI_in[countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_out[FR1].push_back(dataR_in[countR]);
		      countR++;
		    }
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT C " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT C " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes1;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      Full_Stop_Do_Not_Argue(MComms[0]);
      Full_Stop_Do_Not_Argue(MComms[1]);
      FF << " Send DD " << FractalRank << "\n";
      countsa_in.assign(FractalNodes1,0);
      How_Many_Things_To_Send_I(MComms[1],countsa_out,countsa_in);
      dataI_in.clear();
      dataR_in.clear();

      total_in=0;
      total_out=0;
      for(int FR1=0;FR1<FractalNodes1;FR1++)
	{
	  total_out+=countsa_out[FR1];
	  total_in+=countsa_in[FR1];
	}
      FF << " TotalsOther 1 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << endl;
      Send_Data_Somewhere_No_Block(MComms[1],countsa_out,countsa_in,
				   integers+2,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);

      FF << " Send EE " << FractalRank << endl;
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes2);
      dataR_out.resize(FractalNodes2);
      countsa_out.assign(FractalNodes2,0);
      countI=0;
      countR=0;
      try
	{
	  for(int FR1=0;FR1<FractalNodes1;FR1++)
	    {
	      for(int c=0;c<countsa_in[FR1];c++)
		{
		  int FR=dataI_in[countI];
		  countI++;
		  int FR2=FR/FractalNodes01;
		  countsa_out[FR2]++;
		  int FRFrom=dataI_in[countI];
		  countI++;
		  dataI_out[FR2].push_back(FRFrom);
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_out[FR2].push_back(dataI_in[countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_out[FR2].push_back(dataR_in[countR]);
		      countR++;
		    }
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT D " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT D " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes2;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      Full_Stop_Do_Not_Argue(MComms[1]);
      Full_Stop_Do_Not_Argue(MComms[2]);
      FF << " Send FF " << FractalRank << "\n";
      countsa_in.assign(FractalNodes2,0);
      How_Many_Things_To_Send_I(MComms[2],countsa_out,countsa_in);
      dataI_in.clear();
      dataR_in.clear();
      total_in=0;
      total_out=0;
      for(int FR2=0;FR2<FractalNodes2;FR2++)
	{
	  total_out+=countsa_out[FR2];
	  total_in+=countsa_in[FR2];
	}
      FF << " TotalsOther 2 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << endl;
      Send_Data_Somewhere_No_Block(MComms[2],countsa_out,countsa_in,
				   integers+1,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);

      FF << " Send GG " << FractalRank << endl;
      dataI_out.clear();
      dataR_out.clear();
      dataI_out.resize(FractalNodes);
      dataR_out.resize(FractalNodes);
      counts_in.assign(FractalNodes,0);
      countI=0;
      countR=0;
      try
	{
	  for(int FR2=0;FR2<FractalNodes2;FR2++)
	    {
	      for(int c=0;c<countsa_in[FR2];c++)
		{
		  int FRFrom=dataI_in[countI];
		  countI++;
		  counts_in[FRFrom]++;
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_out[FRFrom].push_back(dataI_in[countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_out[FRFrom].push_back(dataR_in[countR]);
		      countR++;
		    }
		}
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT E " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT E " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      dataI_in.clear();
      dataR_in.clear();
      FF << " Send HH " << FractalRank <<  "\n";
      how_manyI=0;
      how_manyR=0;
      try
	{
	  for(int FR=0;FR<FractalNodes;FR++)
	    {
	      countI=0;
	      countR=0;
	      for(int c=0;c<counts_in[FR];c++)
		{
		  for(int nI=0;nI<integers;nI++)
		    {
		      dataI_in.push_back(dataI_out[FR][countI]);
		      countI++;
		    }
		  for(int nR=0;nR<doubles;nR++)
		    {
		      dataR_in.push_back(dataR_out[FR][countR]);
		      countR++;
		    }
		}
	      how_manyI+=countI;
	      how_manyR+=countR;
	    }
	}
      catch(bad_alloc& ba)
	{
	  cout << " DUMP IT F " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  FF << " DUMP IT F " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	  for(int FR=0;FR<FractalNodes;FR++)
	    FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	  FF << endl;
	  assert(0);
	}
      Full_Stop_Do_Not_Argue();
      FF << " Send II " << FractalRank << " " << how_manyI << " " << how_manyR << "\n";
    }
    //
    void Send_Data_Hypre_Directions(vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
				    vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
				    vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR)
    {
      double aNodes=HypreNodes;
      int HypreNodes0=pow(aNodes-0.5,1.0/3.0)+1.0;
      double a12=HypreNodes/HypreNodes0;
      int HypreNodes1=sqrt(a12-0.5)+1.0;
      int HypreNodes01=HypreNodes0*HypreNodes1;
      int HypreNodes2=HypreNodes/HypreNodes01;
      int HypreRank0=HypreRank % HypreNodes0;
      int HypreRank1=(HypreRank/HypreNodes0) % HypreNodes1;
      int HypreRank2=HypreRank/HypreNodes01;
      //      int HypreNodesBox=HypreNodes01*HypreNodes2;
      //      int extras=HypreNodes-HypreNodesBox;
      //      int ExtraLines=extras/HypreNodes0;
      //      int ExtraNodes=extras % HypreNodes0;
      int HypreLong2=how_many_nodes(HComms[2]);
      int HypreLong1=how_many_nodes(HComms[1]);
      int HypreLong0=how_many_nodes(HComms[0]);
      if(HypreRank == 0)
	{
	  cout << " AHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;
	  cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
	}
      vector < vector <int> > dataIa_out(HypreLong2);
      vector < vector <double> > dataRa_out(HypreLong2);
      vector <int>dataIa_in;
      vector <double>dataRa_in;
      vector <int>countsa_out(HypreLong2,0);
      vector <int>countsa_in(HypreLong2);
      int H2start=0;
      for(int HR=0;HR<HypreNodes;HR++)
	{
	  int HR2=HR/HypreNodes01;
	  if(HR2 == HypreNodes2)
	    {
	      H2start=(H2start+1) % HypreNodes2;
	      HR2=H2start;
	    }
	  countsa_out[HR2]+=counts_out[HR];
	  int nIdata=0;
	  int nRdata=0;
	  for(int ni=0;ni<counts_out[HR];ni++)
	    {
	      dataIa_out[HR2].push_back(HR);
	      for(int ints=0;ints<integers;ints++)
		{
		  dataIa_out[HR2].push_back(dataI_out[HR][nIdata]);
		  nIdata++;
		}
	      for(int reals=0;reals<doubles;reals++)
		{
		  dataRa_out[HR2].push_back(dataR_out[HR][nRdata]);
		  nRdata++;
		}
	    }
	}
      //      cout << " BHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataI_out.clear();
      dataR_out.clear();
      How_Many_Things_To_Send_I(HComms[2],countsa_out,countsa_in);
      //      cout << " CHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";

      dataIa_in.clear();
      dataRa_in.clear();

      Send_Data_Somewhere_No_Block(HComms[2],countsa_out,countsa_in,
				   integers+1,doubles,
				   dataIa_out,dataIa_in,how_manyI,
				   dataRa_out,dataRa_in,how_manyR);
      //      cout << " DHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataIa_out.clear();
      dataRa_out.clear();
      dataIa_out.resize(HypreLong1);
      dataRa_out.resize(HypreLong1);
      countsa_out.assign(HypreLong1,0);
      int countI=0;
      int countR=0;
      for(int HR2=0;HR2<HypreLong2;HR2++)
	{
	  int HRFrom=HypreRank+(HR2-HypreRank2)*HypreNodes01;
	  for(int c=0;c<countsa_in[HR2];c++)
	    {
	      int HR=dataIa_in[countI];
	      countI++;
	      int HR1=(HR/HypreNodes0) % HypreNodes1;
	      countsa_out[HR1]++;
	      dataIa_out[HR1].push_back(HR);
	      dataIa_out[HR1].push_back(HRFrom);
	      for(int nI=0;nI<integers;nI++)
		{
		  dataIa_out[HR1].push_back(dataIa_in[countI]);
		  countI++;
		}
	      for(int nR=0;nR<doubles;nR++)
		{
		  dataRa_out[HR1].push_back(dataRa_in[countR]);
		  countR++;
		}
	    }
	}
      Full_Stop_Do_Not_Argue(HComms[2]);
      Full_Stop_Do_Not_Argue(HComms[1]);
      //      cout << " EHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      countsa_in.assign(HypreLong1,0);
      How_Many_Things_To_Send_I(HComms[1],countsa_out,countsa_in);
      dataIa_in.clear();
      dataRa_in.clear();
      //      cout << " FHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";

      Send_Data_Somewhere_No_Block(HComms[1],countsa_out,countsa_in,
				   integers+2,doubles,
				   dataIa_out,dataIa_in,how_manyI,
				   dataRa_out,dataRa_in,how_manyR);
      //      cout << " GHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";

      dataIa_out.clear();
      dataRa_out.clear();
      dataIa_out.resize(HypreLong0);
      dataRa_out.resize(HypreLong0);
      countsa_out.assign(HypreLong0,0);
      countI=0;
      countR=0;
      for(int HR1=0;HR1<HypreLong1;HR1++)
	{
	  for(int c=0;c<countsa_in[HR1];c++)
	    {
	      int HR=dataIa_in[countI];
	      countI++;
	      int HR0=HR % HypreNodes0;
	      countsa_out[HR0]++;
	      int HRFrom=dataIa_in[countI];
	      countI++;
	      dataIa_out[HR0].push_back(HR);
	      dataIa_out[HR0].push_back(HRFrom);
	      for(int nI=0;nI<integers;nI++)
		{
		  dataIa_out[HR0].push_back(dataIa_in[countI]);
		  countI++;
		}
	      for(int nR=0;nR<doubles;nR++)
		{
		  dataRa_out[HR0].push_back(dataRa_in[countR]);
		  countR++;
		}
	    }
	}
      //      cout << " HHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      Full_Stop_Do_Not_Argue(HComms[1]);
      Full_Stop_Do_Not_Argue(HComms[0]);
      countsa_in.assign(HypreLong0,0);
      How_Many_Things_To_Send_I(HComms[0],countsa_out,countsa_in);
      //      cout << " IHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataIa_in.clear();
      dataRa_in.clear();
      Send_Data_Somewhere_No_Block(HComms[0],countsa_out,countsa_in,
				   integers+2,doubles,
				   dataIa_out,dataIa_in,how_manyI,
				   dataRa_out,dataRa_in,how_manyR);
      //      cout << " JHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataIa_out.clear();
      dataRa_out.clear();
      dataIa_out.resize(HypreLong2);
      dataRa_out.resize(HypreLong2);
      countsa_out.assign(HypreLong2,0);
      countI=0;
      countR=0;
      for(int HR0=0;HR0<HypreLong0;HR0++)
	{
	  for(int c=0;c<countsa_in[HR0];c++)
	    {
	      int HR=dataIa_in[countI];
	      countI++;
	      int HR2=HR/HypreNodes01;
	      countsa_out[HR2]++;
	      int HRFrom=dataIa_in[countI];
	      countI++;
	      dataIa_out[HR2].push_back(HRFrom);
	      for(int nI=0;nI<integers;nI++)
		{
		  dataIa_out[HR2].push_back(dataIa_in[countI]);
		  countI++;
		}
	      for(int nR=0;nR<doubles;nR++)
		{
		  dataRa_out[HR2].push_back(dataRa_in[countR]);
		  countR++;
		}
	    }
	}
      //      cout << " KHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      Full_Stop_Do_Not_Argue(HComms[0]);
      Full_Stop_Do_Not_Argue(HComms[2]);
      countsa_in.assign(HypreLong2,0);
      How_Many_Things_To_Send_I(HComms[2],countsa_out,countsa_in);
      //      cout << " LHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataIa_in.clear();
      dataRa_in.clear();
      Send_Data_Somewhere_No_Block(HComms[2],countsa_out,countsa_in,
				   integers+1,doubles,
				   dataIa_out,dataIa_in,how_manyI,
				   dataRa_out,dataRa_in,how_manyR);
      //      cout << " MHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataIa_out.clear();
      dataRa_out.clear();
      dataIa_out.resize(HypreNodes);
      dataRa_out.resize(HypreNodes);
      counts_in.assign(HypreNodes,0);
      countI=0;
      countR=0;
      for(int HR2=0;HR2<HypreLong2;HR2++)
	{
	  for(int c=0;c<countsa_in[HR2];c++)
	    {
	      int HRFrom=dataIa_in[countI];
	      countI++;
	      counts_in[HRFrom]++;
	      for(int nI=0;nI<integers;nI++)
		{
		  dataIa_out[HRFrom].push_back(dataIa_in[countI]);
		  countI++;
		}
	      for(int nR=0;nR<doubles;nR++)
		{
		  dataRa_out[HRFrom].push_back(dataRa_in[countR]);
		  countR++;
		}
	    }
	}
      //      cout << " KHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
      dataI_in.clear();
      dataR_in.clear();
      how_manyI=0;
      how_manyR=0;
      for(int HR=0;HR<HypreNodes;HR++)
	{
	  countI=0;
	  countR=0;
	  for(int c=0;c<counts_in[HR];c++)
	    {
	      for(int nI=0;nI<integers;nI++)
		{
		  dataI_in.push_back(dataIa_out[HR][countI]);
		  countI++;
		}
	      for(int nR=0;nR<doubles;nR++)
		{
		  dataR_in.push_back(dataRa_out[HR][countR]);
		  countR++;
		}
	    }
	  how_manyI+=countI;
	  how_manyR+=countR;
	}
      Full_Stop_Do_Not_Argue(HypreWorld);
      //      cout << " LHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
      //      cout << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    }
    void MPI_MYTest(int which,int test)
    {
      if(test == MPI_SUCCESS)
      	return;
      fprintf(p_file->PFFractalMemory," MPI Error %d %d %d %d %d %d %d %d \n",which,test,
	      MPI_ERR_COMM,MPI_ERR_TYPE,MPI_ERR_COUNT,MPI_ERR_TAG,MPI_ERR_RANK,MPI_SUCCESS);
    }
    template <class T> void my_AllgatherI(vector <T>& paramsend,vector <T>& paramrecv,int nsend)
    {
      int ROOT=FractalNodes/2;
      if(sizeof(T) == sizeof(int))
	{
	  MPI_Gather(&(*(paramsend.begin())),nsend,MPI_INT,&(*(paramrecv.begin())),nsend,MPI_INT,ROOT,FractalWorld);
	  MPI_Bcast(&(*paramrecv.begin()),nsend*FractalNodes,MPI_INT,ROOT,FractalWorld);
	}
      else
	{
	  MPI_Gather(&(*(paramsend.begin())),nsend,MPI_LONG,&(*(paramrecv.begin())),nsend,MPI_LONG,ROOT,FractalWorld);
	  MPI_Bcast(&(*paramrecv.begin()),nsend*FractalNodes,MPI_LONG,ROOT,FractalWorld);
	}
    }
    void my_AllgatherR(vector <double>& paramsend,vector <double>& paramrecv,int nsend)
    {
      int ROOT=FractalNodes/2;
      MPI_Gather(&(*(paramsend.begin())),nsend,MPI_DOUBLE,&(*(paramrecv.begin())),nsend,MPI_DOUBLE,ROOT,FractalWorld);
      Send_DOUBLE_from_ROOT(paramrecv,nsend*FractalNodes,ROOT);
    }
    void calc_total_particles(const int& NP)
    {
      vector <long int> particles(1);
      particles[0]=NP;
      Find_Sum_LONG_INT(particles,1);
      number_particles_total=particles[0];
    }
    void Find_Max_INT(vector <int>& integers,int how_long)
    {
      int ROOT=FractalNodes/2;
      Find_Max_INT_to_ROOT(integers,how_long,ROOT);
      Send_INT_from_ROOT(integers,how_long,ROOT);
    }
    void Find_Max_DOUBLE(vector <double>& doubles,int how_long)
    {
      int ROOT=FractalNodes/2;
      Find_Max_DOUBLE_to_ROOT(doubles,how_long,ROOT);
      Send_DOUBLE_from_ROOT(doubles,how_long,ROOT);
    }
    void Find_Sum_LONG_INT(vector <long int>& integers,int how_long)
    {
      int ROOT=FractalNodes/2;
      Find_Sum_LONG_INT_to_ROOT(integers,how_long,ROOT);
      Send_LONG_INT_from_ROOT(integers,how_long,ROOT);
    }
    void Find_Sum_INT(vector <int>& integers,int how_long)
    {
      int ROOT=FractalNodes/2;
      Find_Sum_INT_to_ROOT(integers,how_long,ROOT);
      Send_INT_from_ROOT(integers,how_long,ROOT);
    }
    void Find_Sum_DOUBLE(vector <double>& doubles,int how_long)
    {
      int ROOT=FractalNodes/2;
      Find_Sum_DOUBLE_to_ROOT(doubles,how_long,ROOT);
      Send_DOUBLE_from_ROOT(doubles,how_long,ROOT);
    }
    void Send_INT_from_ROOT(int* numbers,const int& how_long,const int& ROOT)
    {
      MPI_Bcast(numbers,how_long,MPI_INT,ROOT,FractalWorld);
    }
    void Find_Max_INT_to_ROOT(vector <int>& numbers,int how_long,int ROOT)
    {
      vector <int> maxi(how_long);
      MPI_Reduce(&(*numbers.begin()),&(*maxi.begin()),how_long,MPI_INT,MPI_MAX,ROOT,FractalWorld);
      numbers=maxi;
    }
    void Find_Max_DOUBLE_to_ROOT(vector <double>& numbers,int how_long,int ROOT)
    {
      vector <double> maxr(how_long);
      MPI_Reduce(&(*numbers.begin()),&(*maxr.begin()),how_long,MPI_DOUBLE,MPI_MAX,ROOT,FractalWorld);
      numbers=maxr;
    }
    void Find_Sum_INT_to_ROOT(vector <int>& numbers,int how_long,int ROOT)
    {
      vector <int> sumup(how_long);
      MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_INT,MPI_SUM,ROOT,FractalWorld);
      numbers=sumup;
    }
    void Find_Sum_LONG_INT_to_ROOT(vector <long int>& numbers,int how_long,int ROOT)
    {
      vector <long int> sumup(how_long);
      MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_LONG,MPI_SUM,ROOT,FractalWorld);
      numbers=sumup;
    }
    void Find_Sum_DOUBLE_to_ROOT(vector <double>& numbers,int how_long,int ROOT)
    {
      vector <double> sumup(how_long);
      MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_DOUBLE,MPI_SUM,ROOT,FractalWorld);
      numbers=sumup;
    }
    void Send_INT_from_ROOT(vector <int>& numbers,int how_long,int ROOT)
    {
      MPI_Bcast(&(*numbers.begin()),how_long,MPI_INT,ROOT,FractalWorld);
    }
    void Send_LONG_INT_from_ROOT(vector <long int>& numbers,int how_long,int ROOT)
    {
      MPI_Bcast(&(*numbers.begin()),how_long,MPI_LONG,ROOT,FractalWorld);
    }
    void Send_DOUBLE_from_ROOT(vector <double>& numbers,int how_long,int ROOT)
    {
      MPI_Bcast(&(*numbers.begin()),how_long,MPI_DOUBLE,ROOT,FractalWorld);
    }
    void Full_Stop()
    {
      if(time_trial)
	MPI_Barrier(FractalWorld);
    }
    void Full_Stop_Do_Not_Argue()
    {
      MPI_Barrier(FractalWorld);
    }
    void Full_Stop(MPI_Comm& World)
    {
      if(time_trial)
	MPI_Barrier(World);
    }
    void Full_Stop_Do_Not_Argue(MPI_Comm& World)
    {
      MPI_Barrier(World);
    }
    void zeroR()
    {
      std::fill(potR,potR+2*total_memory,0.0);
    }
    void zeroR(double grail)
    {
      std::fill(potR,potR+2*total_memory,grail);
    }
    double Clock()
    {
      return MPI_Wtime();
    }
    int MyHypreRank()
    {
      int rank;
      if(HypreWorld == MPI_COMM_NULL)
	return -1;
      MPI_Comm_rank(HypreWorld,&rank);
      return rank;      
    }
    void HypreGroupCreate(vector <int>& ranks)
    {
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Group_incl(FractalGroup, HypreNodes, &(*(ranks.begin())), &HypreGroup);
      MPI_Comm_create(FractalWorld, HypreGroup, &HypreWorld);
      HypreRank=MyHypreRank();
      if(!IAmAHypreNode)
	return;
      make_MPI_Hypre_Groups();
    }
    void HypreGroupFree()
    {
      if(!IAmAHypreNode)
      	return;
      MPI_Group_free(&HypreGroup);
      MPI_Comm_free(&HypreWorld);
      HypreWorld=MPI_COMM_NULL;
      HypreGroups3Free();
    }
    void HypreGroups3Free()
    {
      if(!IAmAHypreNode)
      	return;
      for(int ni=0;ni<3;ni++)
	{
	  MPI_Group_free(&HG[ni]);
	  MPI_Comm_free(&HComms[ni]);
	  HComms[ni]=MPI_COMM_NULL;
	}
    }
    void createFractalWorld(MPI_Comm& World,vector <int>& dims)
    {
      MPI_Group_free(&FractalGroup);
      MPI_Comm_free(&FractalWorld);
      FractalNodes=how_many_nodes(World);
      int* Ranks=new int[FractalNodes];
      for(int ni=0;ni<FractalNodes;ni++)
	Ranks[ni]=ni;
      MPI_Group WorldGroup;
      MPI_Comm_group(World,&WorldGroup);
      MPI_Group_incl(WorldGroup, FractalNodes, Ranks, &FractalGroup);
      MPI_Comm_create(World, FractalGroup, &FractalWorld);
      delete [] Ranks;
      dims[0]=max(dims[0],0);
      dims[1]=max(dims[1],0);
      dims[2]=max(dims[2],0);
      MPI_Dims_create(FractalNodes,3,&(*dims.begin()));
      FractalNodes0=dims[0];
      FractalNodes1=dims[1];
      FractalNodes2=dims[2];
    }
  };
}
#endif
