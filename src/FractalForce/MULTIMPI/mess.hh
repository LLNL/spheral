#ifndef _Mess_Defined_
#define _Mess_Defined_
namespace FractalSpace
{
  typedef ptrdiff_t pint;
  class Mess{
  public:
    int FractalRank;
    int FractalNodes;
    int FFTRank;
    int FFTNodes;
    int HypreRank;
    int HypreNodes;
    int number_particles_total;
    Particle* parts_tmp;
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
    vector <int>return_Slice_pos;
    vector <int>return_group;
    vector <int>return_point;
    vector <int>return_node;
    vector <int>what_Slice_point;
    vector <double> green;
    MPI_Comm FractalWorld;
    MPI_Comm FFTWorld;
    MPI_Comm HypreWorld;
    MPI_Group FractalGroup;
    MPI_Group FFTGroup;
    MPI_Group HypreGroup;
    bool IAmAnFFTNode;
    bool IAmAHypreNode;
    bool time_trial;
    bool standalone;
    File* p_file;
    double WallTime;
    Mess():
      FractalRank(0),
      FractalNodes(1),
      FFTRank(0),
      FFTNodes(1234567),
      HypreRank(0),
      HypreNodes(0),
      number_particles_total(-1),
      start_x(0),
      length_x(-1),
      total_memory(-1),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      standalone(true)
    {
      WallTime=Clock();
      cout << " Empty Mess " << endl;
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP,const int& FN,MPI_Comm& FW):
      FractalRank(0),
      FractalNodes(1),
      FFTNodes(FN),
      HypreRank(0),
      HypreNodes(0),
      number_particles_total(-1),
      start_x(0),
      length_x(GR),
      total_memory(1),
      FractalWorld(FW),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true)
    {
      cout << " Making a Mess with parameters" << endl;
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
      cout << " starting to clean up a mess " << FractalRank << endl;
      FFTWFinal();
      if(standalone)
	MPIFinal();
      cout << " cleaned up a mess " << FractalRank << endl;
    }
    void MPIStartup()
    {
      cout << " Into MPIStartup " << endl;
      int knights;
      MPI_Initialized(&knights);
      if(!knights)
	MPI_Init(NULL,NULL);
      FractalWorld=MPI_COMM_WORLD;
      FFTWorld=FractalWorld;
      HypreWorld=FractalWorld;
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Comm_group(FFTWorld,&FFTGroup);
      //      MPI_Comm_group(HypreWorld,&HypreGroup);
      FractalRank=what_is_my_rank(); 
      FractalNodes=how_many_nodes();
      FFTRank=FractalRank;
      FFTNodes=min(FFTNodes,FractalNodes);
      HypreRank=FractalRank;
      HypreNodes=FractalNodes;
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
    int what_is_my_FFT_rank()
    {
      int rank;
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
      cout << " FFTWStartup a " << FractalRank << " " << FFTRank << " " << FFTNodes << endl;
      doFFTWorld(length_1,periodic);
      fftw_mpi_init();
      cout << " FFTWStartup b " << FractalRank << " " << FFTRank << " " << FFTNodes << endl;
      if(!IAmAnFFTNode) return;
      if(periodic)
	{
	  const pint Length_c=(Length_1+2)/2;
	  total_memory=fftw_mpi_local_size_3d(Length_1,Length_1,Length_c,FFTWorld,&length_x,&start_x);
	  create_potRC();
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(Length_1,Length_1,Length_1,potR,potC,FFTWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(Length_1,Length_1,Length_1,potC,potR,FFTWorld,FFTW_MEASURE);
	}
      else
	{
	  const pint Length_11=Length_1+1;
	  const pint Length_2=2*Length_1;
	  const double g_c=pow(static_cast<double>(Length_1),-5)/8.0;
	  cout << " g_c= " << g_c << " " << FractalRank << endl;
	  total_memory=fftw_mpi_local_size_3d(Length_2,Length_2,Length_11,FFTWorld,&length_x,&start_x);
	  cout << " total_memory " << FractalRank << " " << total_memory << " " << length_x << " " << start_x << endl;
	  cout << "mess haha " << FractalRank << " " << green.max_size() << endl;
	  green.resize(length_x*Length_11*Length_11);
	  cout << " mess a " << FractalRank << endl;
	  create_potRC();
	  cout << " mess b " << FractalRank << endl;
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(Length_2,Length_2,Length_2,potR,potC,FFTWorld,FFTW_MEASURE);
	  cout << " mess c " << FractalRank << endl;
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(Length_2,Length_2,Length_2,potC,potR,FFTWorld,FFTW_MEASURE);
	  cout << " mess d " << FractalRank << endl;
	  zeroR();
	  cout << " mess e " << FractalRank << endl;
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
      cout << " messya " << FractalRank << " " << length_1 << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << IAmAnFFTNode << endl;
      FFTRank=FractalRank;
      int maxFFT=length_1/2;
      if(!periodic)
	maxFFT*=2;
      FFTNodes=min(FFTNodes,FractalNodes);
      FFTNodes=min(FFTNodes,maxFFT);
      FFTNodes=(FFTNodes/2)*2;
      IAmAnFFTNode=FFTRank < FFTNodes;
      cout << " messyb " << FractalRank << " " << length_1 << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << maxFFT << " " << IAmAnFFTNode << endl;
      int* ranks=new int[FFTNodes];
      for(int ni=0;ni < FFTNodes;ni++)
	ranks[ni]=ni;
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Group_incl(FractalGroup, FFTNodes, ranks, &FFTGroup);
      MPI_Comm_create(FractalWorld, FFTGroup, &FFTWorld);
      delete [] ranks;
      ranks=0;
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
      cout << " messyc " << FractalRank << " " << length_1 << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << maxFFT << " " << IAmAnFFTNode << endl;
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
    void dumpR(ofstream& FILE,const int& length,const bool& test)
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
      int paramsend[2]={(int)start_x,(int)(start_x+length_x-1)};
      //      int* paramrecv=(int*)malloc(2*FractalNodes*sizeof(int));
      int* paramrecv=new int[2*FractalNodes];
      int length_1=length_a;
      if(!periodic)
	length_1=2*length_a;
      paramsend[0]=start_x;
      paramsend[1]=start_x+length_x-1;
      cout << "calc_fftwa " << FFTRank << " " << start_x << " " << length_x << endl;
      //      if(IAmAnFFTNode)
      MPI_Allgather(paramsend,2,MPI_INT,paramrecv,2,MPI_INT,FractalWorld);
      //      MPI_Barrier(FractalWorld);
      cout << "calc_fftwb " << FFTRank << " " << start_x << " " << length_x << endl;
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
	  cout << " slices " << FFTRank << " " << Slices[FR][0] << " " << Slices[FR][1] << " " << FR << " " << FractalRank << endl;
	}
      delete [] paramrecv;
      //      free(paramrecv);
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
		cout << " success " << FractalRank << " " << FFTRank << " " << nx << " " << WhichSlice[nx] << endl;
	    }
	}
      for(int ni=0;ni<length_1;ni++)
	cout << "whichslice " << FFTRank << " " << ni << " " << WhichSlice[ni] << endl;
      assert(allok);
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
      assert(p_file);
      cout << " into how many " << FractalNodes << " " << p_file << endl;
      //      ofstream& FF=p_file->FileFractal;
      int* counts_out=new int[FractalNodes];
      int* counts_in=new int[FractalNodes];
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  counts_out[FR]=counts_out_send[FR];
	  //	  FF << "sending a " << FractalRank << " " << FR << " " << counts_out[FR] << endl;
	}
      for(int FR=0;FR<FractalNodes;FR++)
	MPI_Gather(&counts_out[FR],1,MPI_INT,counts_in,1,MPI_INT,FR,FractalWorld);
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  counts_in_send[FR]=counts_in[FR];
	  //	  FF << "sending b " << FractalRank << " " << FR << " " << counts_out[FR] << " " << counts_in[FR] << endl;
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
      //      cout << " howmanyIR " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
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
      //      cout << " how many " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
      if(integers > 0)
	{
	  dataI_in_send.resize(how_manyI);
	  //      std::copy(dataI_in,dataI_in+how_manyI,dataI_in_send.begin());
	  for(int ni=0;ni<how_manyI;ni++)
	    dataI_in_send[ni]=dataI_in[ni];
	}
      if(doubles > 0)
	{
	  dataR_in_send.resize(how_manyR);
	  //      std::copy(dataR_in,dataR_in+how_manyR,dataR_in_send.begin());
	  for(int ni=0;ni<how_manyR;ni++)
	    dataR_in_send[ni]=dataR_in[ni];
	}
      delete [] dataI_in;
      delete [] dataR_in;
      dataI_in=0;
      dataR_in=0;
    }
    void Send_Data_Somewhere_Faster(vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
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
	  if(FR != FractalRank)
	    maxcount=max(maxcount,counts_out_send[FR]);
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
      countsI_out[FractalRank]=0;
      countsR_out[FractalRank]=0;
      int* DataI_out=new int[max(maxcount*integers,1)];
      double* DataR_out=new double[max(maxcount*doubles,1)];
      int* dataI_in=new int[how_manyI];
      double* dataR_in=new double[how_manyR];
      //      cout << " howmanyIR " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
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
      //      cout << " how many " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
      if(integers > 0)
	{
	  dataI_in_send.resize(how_manyI);
	  //      std::copy(dataI_out[FractalRank].begin(),dataI_out[FractalRank].begin()+extraI,dataI_in.begin()+startI);
	  for(int ni=0;ni<extraI;ni++)
	    dataI_in[ni+startI]=dataI_out[FractalRank][ni];
	  //      std::copy(dataI_in,dataI_in+how_manyI,dataI_in_send.begin());
	  for(int ni=0;ni<how_manyI;ni++)
	    dataI_in_send[ni]=dataI_in[ni];
	}
      if(doubles > 0)
	{
	  dataR_in_send.resize(how_manyR);
	  //      std::copy(dataR_out[FractalRank].begin(),dataR_out[FractalRank].begin()+extraR,dataR_in.begin()+startR);
	  for(int ni=0;ni<extraR;ni++)
	    dataR_in[ni+startR]=dataR_out[FractalRank][ni];
	  //      std::copy(dataR_in,dataR_in+how_manyR,dataR_in_send.begin());
	  for(int ni=0;ni<how_manyR;ni++)
	    dataR_in_send[ni]=dataR_in[ni];
	}
      delete [] dataI_in;
      delete [] dataR_in;
      dataI_in=0;
      dataR_in=0;
    }
    void Send_Data_Somewhere_No_Block(vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
			     vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
			     vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
    {
      MPI_Barrier(FractalWorld);
      ofstream& FF=p_file->FileFractal;
      FF << " in messy faster " << endl;
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
	  FF << "outs " << FR << " " << countsI_in[FR] << " " << countsI_out[FR] << " " << countsR_in[FR] << " " << countsR_out[FR] << " ";
	  FF << displsI[FR] << " " << displsR[FR] << endl;
	}
      how_manyI=displsI[FractalNodes-1]+counts_in_send[FractalNodes-1]*integers;
      how_manyR=displsR[FractalNodes-1]+counts_in_send[FractalNodes-1]*doubles;
      int* DataI_out=new int[max(maxcount*integers,1)];
      double* DataR_out=new double[max(maxcount*doubles,1)];
      int* dataI_in=new int[how_manyI];
      double* dataR_in=new double[how_manyR];
      int tagI=0;
      int tagR=1;
      MPI_Request* requestI=new MPI_Request[FractalNodes];
      MPI_Request* requestR=new MPI_Request[FractalNodes];
      MPI_Status* statusI=new MPI_Status[FractalNodes];
      MPI_Status* statusR=new MPI_Status[FractalNodes];
      FF << " howmanyIR " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  if(FR == FractalRank)
	    {
	      int disp=displsI[FR];
	      for(int ni=0;ni<countsI_out[FR];ni++)
		dataI_in[disp+ni]=dataI_out[FR][ni];
	      disp=displsR[FR];
	      for(int ni=0;ni<countsR_out[FR];ni++)
		dataR_in[disp+ni]=dataR_out[FR][ni];
	      continue;
	    }
	  if(integers > 0)
	    {
	      for(int ni=0;ni<countsI_out[FR];ni++)
		DataI_out[ni]=dataI_out[FR][ni];
	      MPI_Isend(DataI_out,countsI_out[FR],MPI_INT,FR,tagI,FractalWorld,&requestI[FR]);
	    }
	  if(doubles > 0)
	    {
	      for(int ni=0;ni<countsR_out[FR];ni++)
		DataR_out[ni]=dataR_out[FR][ni];
	      MPI_Isend(DataR_out,countsR_out[FR],MPI_DOUBLE,FR,tagR,FractalWorld,&requestR[FR]);
	    }
	}
      //      int flagI=-1;
      //      int flagR=-1;
      int countI=-1;
      int countR=-1;
      MPI_Status statPI;
      MPI_Status statPR;
      for(int FRin=0;FRin<FractalNodes-1;FRin++)
	{
	  if(integers > 0)
	    {
	      MPI_Probe(MPI_ANY_SOURCE,tagI,FractalWorld,&statPI);
	      //	      MPI_Iprobe(MPI_ANY_SOURCE,tagI,FractalWorld,&flagI,&statPI);
	      //	      assert(flagI);
	      MPI_Get_count(&statPI,MPI_INT,&countI);
	      int FR=statPI.MPI_SOURCE;
	      assert(countI==countsI_in[FR]);
	      MPI_Recv(&dataI_in[displsI[FR]],countI,MPI_INT,FR,tagI,FractalWorld,&statusI[FR]);
	    }
	  if(doubles > 0)
	    {
	      MPI_Probe(MPI_ANY_SOURCE,tagR,FractalWorld,&statPR);
	      //	      MPI_Iprobe(MPI_ANY_SOURCE,tagR,FractalWorld,&flagR,&statPR);
	      //	      assert(flagR);
	      MPI_Get_count(&statPR,MPI_DOUBLE,&countR);
	      int FR=statPR.MPI_SOURCE;
	      assert(countR==countsR_in[FR]);
	      MPI_Recv(&dataR_in[displsR[FR]],countR,MPI_DOUBLE,FR,tagR,FractalWorld,&statusR[FR]);
	    }
	}
      if(integers > 0)
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    if(FR != FractalRank)
	      MPI_Wait(&requestI[FR],&statusI[FR]);
	  }
      if(doubles > 0)
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    if(FR != FractalRank)
	      MPI_Wait(&requestR[FR],&statusR[FR]);
	  }
      delete [] requestI;
      delete [] requestR;
      delete [] statusI;
      delete [] statusR;
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
      //      cout << " how many " << FractalRank << " " << how_manyI << " " << how_manyR << endl;
      if(integers > 0)
	{
	  dataI_in_send.resize(how_manyI);
	  //      std::copy(dataI_in,dataI_in+how_manyI,dataI_in_send.begin());
	  for(int ni=0;ni<how_manyI;ni++)
	    dataI_in_send[ni]=dataI_in[ni];
	}
      if(doubles > 0)
	{
	  dataR_in_send.resize(how_manyR);
	  //      std::copy(dataR_in,dataR_in+how_manyR,dataR_in_send.begin());
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
      //      std::copy(maxy,maxy+how_long,integers);
      for(int ni=0;ni<how_long;ni++)
	integers[ni]=maxy[ni];
      delete [] maxy;
    }
    void Find_Max_DOUBLE(double* doubles,const int& how_long)
    {
      double* maxy=new double[how_long];
      MPI_Allreduce(doubles,maxy,how_long,MPI_DOUBLE,MPI_MAX,FractalWorld);
      //      std::copy(maxy,maxy+how_long,doubles);
      for(int ni=0;ni<how_long;ni++)
	doubles[ni]=maxy[ni];
      delete [] maxy;
    }
    void Find_Sum_INT(int* integers,const int& how_long)
    {
      int* sumup=new int[how_long];
      MPI_Allreduce(integers,sumup,how_long,MPI_INT,MPI_SUM,FractalWorld);
      //      std::copy(sumup,sumup+how_long,integers);
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
      //      std::copy(doublesb,doublesb+how_long,doubles.begin());
      for(int ni=0;ni<how_long;ni++)
	doubles[ni]=doublesb[ni];
      delete [] doublesa;
      delete [] doublesb;
    }
    void Find_Sum_INT_to_ROOT(int* numbers,const int& how_long,const int& ROOT)
    {
      int* sumup=new int[how_long];
      MPI_Reduce(numbers,sumup,how_long,MPI_INT,MPI_SUM,ROOT,FractalWorld);
      //      std::copy(sumup,sumup+how_long,numbers);
      for(int ni=0;ni<how_long;ni++)
	numbers[ni]=sumup[ni];
      delete [] sumup;
    }
    void Send_INT_from_ROOT(int* numbers,const int& how_long,const int& ROOT)
    {
      MPI_Bcast(numbers,how_long,MPI_INT,ROOT,FractalWorld);
    }
    void Full_Stop()
    {
      if(time_trial)
	MPI_Barrier(FractalWorld);
    }
    void zeroR()
    {
      std::fill(potR,potR+2*total_memory,0.0);
      //      for(pint ni=0;ni<2*total_memory;ni++)
      //	potR[ni]=0.0;
    }
    void zeroR(double grail)
    {
      std::fill(potR,potR+2*total_memory,grail);
      //      for(pint ni=0;ni<2*total_memory;ni++)
      //	potR[ni]=grail;
    }
    double Clock()
    {
      return MPI_Wtime();
      //      return clock();
    }
    void HypreGroupCreate(vector <int>& ranks)
    {
      int* Ranks=new int[HypreNodes];
      for(int ni=0;ni<HypreNodes;ni++)
	Ranks[ni]=ranks[ni];
      MPI_Comm_group(FractalWorld,&FractalGroup);
      MPI_Group_incl(FractalGroup, HypreNodes, Ranks, &HypreGroup);
      MPI_Comm_create(FractalWorld, HypreGroup, &HypreWorld);
      delete [] Ranks;
    }
    void HypreFree()
    {
      if(!IAmAHypreNode)
	return;
      MPI_Group_free(&HypreGroup);
      MPI_Comm_free(&HypreWorld);
    }
  };
}
#endif
