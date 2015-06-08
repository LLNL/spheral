#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  bool Mess::IAMROOT;
  void Mess::MPIStartup()
  {
    //      cerr << " Into MPIStartup " << "\n";
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
    //    if(FractalRank == 0)
    //      cerr << " initialized MPI " << FractalRank << " " << FractalNodes << "\n";
  }
  void Mess::MPIStartup(const bool& PR,int& FR0,int& FR1,int& FR2)
  {
    //      int ranky;
    //      MPI_Comm_rank(FractalWorld,&ranky);
    //      cerr << " Into MPIStartup A " << ranky << "\n";
    int knights;
    MPI_Initialized(&knights);
    if(!knights)
      MPI_Init(NULL,NULL);
    //      cerr << " Into MPIStartup B " << ranky << "\n";
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
    //      cerr << " initialized MPI " << FractalRank << " " << FractalNodes << "\n";
  }
  void Mess::MPIFinal() const
  {
    int knights;
    MPI_Finalized(&knights);
    if(!knights)
      MPI_Finalize();
  }
  int Mess::what_is_my_rank() const
  {
    int rank;
    MPI_Comm_rank(FractalWorld,&rank);
    return rank;
  }
  int Mess::what_is_my_rank(MPI_Comm& World) const
  {
    int rank;
    MPI_Comm_rank(World,&rank);
    return rank;
  }
  int Mess::how_many_nodes(MPI_Comm& World) const
  {
    int size;
    MPI_Comm_size(World,&size);
    return size;
  }
  int Mess::how_many_nodes() const
  {
    int size;
    MPI_Comm_size(FractalWorld,&size);
    return size;
  }
  int Mess::what_is_my_FFT_rank() const
  {
    int rank;
    if(FFTWorld == MPI_COMM_NULL)
      return -1;
    MPI_Comm_rank(FFTWorld,&rank);
    return rank;
  }
  int Mess::how_many_FFT_nodes() const
  {
    int size;
    MPI_Comm_size(FFTWorld,&size);
    return size;
  }
  int Mess::what_is_my_Hypre_rank() const
  {
    if(!IAmAHypreNode)
      return -1;
    int rank;
    MPI_Comm_rank(HypreWorld,&rank);
    return rank;
  }
  int Mess::how_many_Hypre_nodes() const
  {
    if(!IAmAHypreNode)
      return -1;
    int size;
    MPI_Comm_size(HypreWorld,&size);
    return size;
  }
  void Mess::FFTWStartup(const int& length_1,const bool& periodic)
  {
    const pint Length_1=length_1;
    doFFTWorld(length_1,periodic);
    fftw_mpi_init();
    if(!IAmAnFFTNode) return;
    if(periodic)
      {
	const pint Length_c=(Length_1+2)/2;
	total_memory=fftw_mpi_local_size_3d(Length_1,Length_1,Length_c,FFTWorld,&length_x,&start_x);
	//	cerr << " total_memory " << FractalRank << " " << FFTRank << " " << total_memory << " " << length_x << " " << start_x << "\n";
	create_potRC();
	plan_rc=fftw_mpi_plan_dft_r2c_3d(Length_1,Length_1,Length_1,potR,potC,FFTWorld,FFTW_ESTIMATE);
	plan_cr=fftw_mpi_plan_dft_c2r_3d(Length_1,Length_1,Length_1,potC,potR,FFTWorld,FFTW_ESTIMATE);
	free_potRC();
      }
    else
      {
	const pint Length_11=Length_1+1;
	const pint Length_2=2*Length_1;
	double g_c=pow(static_cast<double>(Length_1),-5)/8.0;
	//	  cerr << " g_c= " << g_c << " " << FractalRank << "\n";
	total_memory=fftw_mpi_local_size_3d(Length_2,Length_2,Length_11,FFTWorld,&length_x,&start_x);
	//	cerr << " total_memory " << FractalRank << " " << FFTRank << " " << total_memory << " " << length_x << " " << start_x << " " << g_c << "\n";
	green.resize(length_x*Length_11*Length_11);
	create_potRC();
	plan_rc=fftw_mpi_plan_dft_r2c_3d(Length_2,Length_2,Length_2,potR,potC,FFTWorld,FFTW_ESTIMATE);
	plan_cr=fftw_mpi_plan_dft_c2r_3d(Length_2,Length_2,Length_2,potC,potR,FFTWorld,FFTW_ESTIMATE);
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
  void Mess::doFFTWorld(int how_long,const bool& periodic)
  {
    if(!periodic)
      how_long*=2;
    int each=2;
    bool keep_trying=true;
    while(keep_trying)
      {
	FFTNodes=how_long/each;
	keep_trying=!(how_long % each == 0 && FFTNodes <= FractalNodes);
	each+=2;
	assert(each < 12);
      }
    IAmAnFFTNode=false;
    Franks.clear();
    ItIsAnFFTNode.assign(FractalNodes,false);
    for(int FN=0;FN<FFTNodes;FN++)
      {
	int FRank=(FN*FractalNodes)/FFTNodes;
	ItIsAnFFTNode[FRank]=true;
	Franks.push_back(FRank);
	//	if(FractalRank == 0)
	//	  cerr << " FFTNODES " << FN << " " << Franks.back() << "\n";
	if(FRank == FractalRank)
	  {
	    IAmAnFFTNode=true;
	    FFTRank=FN;
	  }
      }
    MPI_Comm_group(FractalWorld,&FractalGroup);
    MPI_Group_incl(FractalGroup, FFTNodes,&(*Franks.begin()), &FFTGroup);
    MPI_Comm_create(FractalWorld, FFTGroup, &FFTWorld);
    if(!IAmAnFFTNode)
      {
	start_x=9876543;
	length_x=0;
	total_memory=1;
      }
    else
      {
	assert(FFTRank == what_is_my_FFT_rank());
	assert(FFTNodes == how_many_FFT_nodes());
      }
    //    if(FractalRank == 0)
    //      cerr << " messyc " << FractalRank << " " << how_long << " " << periodic << " " << FFTRank << " " << FFTNodes << " " << IAmAnFFTNode << "\n";
  }
  void Mess::FFTWFinal()
  {
    if(IAmAnFFTNode)
      {
	fftw_destroy_plan(plan_rc);
	fftw_destroy_plan(plan_cr);
      }
    fftw_mpi_cleanup();
  }
  void Mess::dumpR(ofstream& FILE,const int& length) const
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
//   void Mess::create_potRC()
//   {
//     size_t sizeR=sizeof(double);
//     size_t sizeC=sizeof(fftw_complex);
//     potR=(double*) fftw_malloc(sizeR*2*total_memory);
//     potC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
//   }
//   void Mess::free_potRC()
//   {
//     fftw_free(potR);
//     fftw_free(potC);
//   }
//   void Mess::create_potR()
//   {
//     size_t sizeR=sizeof(double);
//     potR=(double*) fftw_malloc(sizeR*2*total_memory);
//   }
//   void Mess::free_potR()
//   {
//     fftw_free(potR);
//   }
//   void Mess::create_potRS()
//   {
//     size_t sizeR=sizeof(double);
//     if(IAmPeriodic)
//       potRS=(double*) fftw_malloc(sizeR*2*total_memory);
//     else
//       potRS=(double*) fftw_malloc(sizeR*(glength+1)*(glength+1)*length_x);
//   }
//   void Mess::free_potRS()
//   {
//     fftw_free(potRS);
//   }
//   void Mess::create_potC()
//   {
//     size_t sizeC=sizeof(fftw_complex);
//     potC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
//   }
//   void Mess::free_potC()
//   {
//     fftw_free(potC);
//   }

  void Mess::create_potRC()
  {
    try
      {
	potR=new double[2*total_memory];
      }
    catch(bad_alloc& ba)
      {
	cerr << " BADD potR 1 " << ba.what() << " " << FractalRank << " " << total_memory << endl;
	assert(0);
      }
    try
      {
	potC=new fftw_complex[total_memory];
      }
    catch(bad_alloc& ba)
      {
	cerr << " BADD potC 1 " << ba.what() << " " << FractalRank << " " << total_memory << endl;
	assert(0);
      }
//     size_t sizeR=sizeof(double);
//     size_t sizeC=sizeof(fftw_complex);
//     potR=(double*) fftw_malloc(sizeR*2*total_memory);
//     potC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
  }
  void Mess::free_potRC()
  {
    delete [] potR;
    delete [] potC;
//     fftw_free(potR);
//     fftw_free(potC);
  }
  void Mess::create_potR()
  {
    try
      {
	potR=new double[2*total_memory];
      }
    catch(bad_alloc& ba)
      {
	cerr << " BADD potR 2 " << ba.what() << " " << FractalRank << " " << total_memory << endl;
	assert(0);
      }
//     size_t sizeR=sizeof(double);
//     potR=(double*) fftw_malloc(sizeR*2*total_memory);
  }
  void Mess::free_potR()
  {
    delete [] potR;
//     fftw_free(potR);
  }
  void Mess::create_potRS()
  {
    //    size_t sizeR=sizeof(double);
    if(IAmPeriodic)
      try
	{
	  potRS=new double[2*total_memory];
	}
      catch(bad_alloc& ba)
	{
	  cerr << " BADD potRS 1 " << ba.what() << " " << FractalRank << " " << total_memory << endl;
	  assert(0);
	}
//       potRS=(double*) fftw_malloc(sizeR*2*total_memory);
    else
      try
	{
	  potRS=new double[(glength+1)*(glength+1)*length_x];
	}
      catch(bad_alloc& ba)
	{
	  cerr << " BADD potRS 1 " << ba.what() << " " << FractalRank << " " << total_memory << endl;
	  assert(0);
	}
//       potRS=(double*) fftw_malloc(sizeR*(glength+1)*(glength+1)*length_x);
  }
  void Mess::free_potRS()
  {
    delete [] potRS;
//     fftw_free(potRS);
  }
  void Mess::create_potC()
  {
    try
      {
	potC=new fftw_complex[total_memory];
      }
    catch(bad_alloc& ba)
      {
	cerr << " BADD potC 2 " << ba.what() << " " << FractalRank << " " << total_memory << endl;
	assert(0);
      }
//     size_t sizeC=sizeof(fftw_complex);
//     potC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
  }
  void Mess::free_potC()
  {
    delete [] potC;
//     fftw_free(potC);
  }


  void Mess::fftw_real_to_complex()
  {
    if(IAmAnFFTNode)
      fftw_mpi_execute_dft_r2c(plan_rc,potR,potC);
  }
  void Mess::fftw_complex_to_real()
  {
    if(IAmAnFFTNode)
      fftw_mpi_execute_dft_c2r(plan_cr,potC,potR);
  }
  int Mess::fftw_where(const int& i,const int& j,const int& k,const int& lb,const int& lc) const
  {
    return k+(j+(i-start_x)*lb)*lc;
  }
  void Mess::calc_fftw_Slices(const int& length_a,const bool& periodic)
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
    //      cerr << "calc_fftwa " << FFTRank << " " << start_x << " " << length_x << "\n";
    my_AllgatherI(paramsend,paramrecv,2);
    //      MPI_Allgather(paramsend,2,MPI_INT,paramrecv,2,MPI_INT,FractalWorld);
    //    if(IAmAnFFTNode)
    //      cerr << "calc_fftwb " << FFTRank << " " << FractalRank << " " << start_x << " " << length_x << "\n";
    //	cerr << "calc_fftwb " << FFTRank << " " << FractalRank << " " << start_x << " " << length_x << "\n";
    Slices.resize(FractalNodes); // this is not an error.
    BoxS.resize(FractalNodes); // it must be dimensioned
    BoxSL.resize(FractalNodes); // this way
    for(int FR=0;FR<FFTNodes;FR++)
      {
	int FFTR=Franks[FR];
	Slices[FFTR].resize(2);
	Slices[FFTR][0]=paramrecv[2*FFTR];
	Slices[FFTR][1]=paramrecv[2*FFTR+1];
	BoxS[FFTR].resize(6);
	BoxS[FFTR][0]=Slices[FFTR][0];
	BoxS[FFTR][1]=Slices[FFTR][1];
	BoxS[FFTR][2]=0;
	BoxS[FFTR][3]=length_1-1;
	BoxS[FFTR][4]=0;
	BoxS[FFTR][5]=length_1-1;
	BoxSL[FFTR].resize(3);
	BoxSL[FFTR][0]=length_x;
	BoxSL[FFTR][1]=length_1;
	BoxSL[FFTR][2]=length_1;
	//	if(FFTRank == 0)
	//	  cerr << " slices " << FFTRank << " " << Slices[FFTR][0] << " " << Slices[FFTR][1] << " " << FR << " " << FractalRank << "\n";
	//	    cerr << " slices " << FFTRank << " " << Slices[FR][0] << " " << Slices[FR][1] << " " << FR << " " << FractalRank << "\n";
      }
    WhichSlice.assign(length_1,-10);
    bool allok=true;
    for(int nx=0;nx<length_1;nx++)
      {
	bool success=false;
	for(int S=0;S<FFTNodes;S++)
	  {
	    int FRS=Franks[S];
	    if(nx >= Slices[FRS][0] && nx <= Slices[FRS][1])
	      {
		WhichSlice[nx]=FRS;
		//		  WhichSlice[nx]=S;
		success=true;
		break;
	      }
	  }
	if(!success)
	  {
	    allok=false;
	    for(int nx=0;nx<length_1;nx++)
	      if(FFTRank == 0) cerr << " success " << FractalRank << " " << FFTRank << " " << nx << " " << WhichSlice[nx] << "\n";
	  }
      }
    //    for(int ni=0;ni<length_1;ni++)
    //      if(FFTRank == 0) cerr << "whichslice " << FFTRank << " " << ni << " " << WhichSlice[ni] << "\n";
    assert(allok);
  }
  void Mess::How_Many_On_Nodes(int count,vector <int>& counts)
  {
    counts.resize(FractalNodes);
    vector <int>counts_out;
    counts_out.push_back(count);
    my_AllgatherI(counts_out,counts,1);
  }
  void Mess::How_Many_On_Nodes(long int count,vector <long int>& counts)
  {
    counts.resize(FractalNodes);
    vector <long int>counts_out;
    counts_out.push_back(count);
    my_AllgatherI(counts_out,counts,1);
  }
  void Mess::MAX_Things_To_Send_Receive_I(vector <int>& counts_out_send,vector <int>& counts_in_send,vector <int>& maxSR)
  {
    ofstream& FF=p_file->DUMPS;
    double time0=Clock();
    How_Many_Things_To_Send_I(counts_out_send,counts_in_send);
    double time1=Clock();
    maxSR.clear();
    maxSR.push_back(std::accumulate(counts_out_send.begin(),counts_out_send.end(),0));
    maxSR.push_back(std::accumulate(counts_in_send.begin(),counts_in_send.end(),0));
    double time2=Clock();
    Find_Max_INT(maxSR,2);
    double time3=Clock();
    FF << " SENDREC " << time1-time0 << " " << time2-time1 << " " << time3-time2 << "\n";
  }
  //
  //
  void Mess::How_Many_Things_To_Send_I(vector <int>& counts_out_send,vector <int>& counts_in_send)
  {
    How_Many_Things_To_Send_I(FractalWorld,counts_out_send,counts_in_send);
  }
  void Mess::How_Many_Things_To_Send_I(MPI_Comm& World,
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
  void Mess::Send_Data_Somewhere_No_Block(vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
					  vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
					  vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
  {
    Send_Data_Somewhere_No_Block(FractalWorld,
				 counts_out_send,counts_in_send,integers,doubles,
				 dataI_out,dataI_in_send,how_manyI,
				 dataR_out,dataR_in_send,how_manyR);
  }
  void Mess::Send_Data_Somewhere_No_Block(MPI_Comm& World,
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
	try
	  {
	    dataI_in_send.resize(how_manyI);
	  }
	catch(bad_alloc& ba)
	  {
	    cerr << "BAAD DATAI_IN " << ba.what() << " " << Rank << " " << Nodes << " " << how_manyI << " " << how_manyR << endl;
	    for(int FR=0;FR<Nodes;FR++)
	      cerr << " COUNTSI " << FR << " " << counts_out_send[FR] << " " << counts_in_send[FR] << "\n";
	    assert(0);
	  }
	//	  for(int ni=0;ni<extraI;ni++)
	//	    dataI_in_send[ni+startI]=dataI_out[Rank][ni];
	std::copy(dataI_out[Rank].begin(),dataI_out[Rank].begin()+extraI,dataI_in_send.begin()+startI);
      }
    if(doubles > 0)
      {
	try
	  {
	    dataR_in_send.resize(how_manyR);
	  }
	catch(bad_alloc& ba)
	  {
	    cerr << "BAAD DATAR_IN " << ba.what() << " " << Rank << " " << Nodes << " " << how_manyI << " " << how_manyR << endl;
	    for(int FR=0;FR<Nodes;FR++)
	      cerr << " COUNTSR " << FR << " " << counts_out_send[FR] << " " << counts_in_send[FR] << "\n";
	    assert(0);
	  }
	//	  for(int ni=0;ni<extraR;ni++)
	//	    dataR_in_send[ni+startR]=dataR_out[Rank][ni];
	std::copy(dataR_out[Rank].begin(),dataR_out[Rank].begin()+extraR,dataR_in_send.begin()+startR);
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
  void Mess::make_MPI_Hypre_Groups()
  {
    //      cerr << " making hypre " << FractalRank << " " << HypreRank << "\n";
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
//     if(HypreRank == 0)
//       {
// 	cerr << " MAKE HGD " << HypreNodes << " " << HypreNodes0 << " " << HypreNodes1 << " " << HypreNodes2 << " " << HypreNodes01;
// 	cerr << " " << HypreRank0 <<  " " << HypreRank1 <<  " " << HypreRank2 << " " << ExtraLines << " " << ExtraNodes << "\n";
//       }
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
  void Mess::make_MPI_Groups()
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
  void Mess::Send_Data_Some_How(int tag,vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
				vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
				vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
  {
    Send_Data_Some_How(tag,FractalWorld,counts_out,counts_in,integers,doubles,
		       dataI_out,dataI_in_send,how_manyI,
		       dataR_out,dataR_in_send,how_manyR);
  }
  void Mess::Send_Data_Some_How(int tag,MPI_Comm& World,
				vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
				vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
				vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR)
  {
    fftwTAG=tag;
    //    int Rank=what_is_my_rank(World);
    int Nodes=how_many_nodes(World);
    bool small=Nodes <= MPI_SWITCH;
    bool foreign=World != FractalWorld && World != HypreWorld;
    possibleDANGER=Nodes >= MPI_MAX_COMMS && (fftwTAG == 0 || fftwTAG == 4 || fftwTAG == 7) && World == FractalWorld;
    DANGERlevel=possibleDANGER ? 1:0;
    //    cerr << " tag= " << tag << " " << Rank << " " << Nodes << " " << small << " " << foreign << "\n";
    //    if(Rank == 0)
    //      cerr << " SOMEWHOW " << FractalRank << " " << Nodes << " " << MPI_SWITCH << " " << small << foreign << " ";
    if(small || foreign)
      {
	//	if(Rank == 0)
	//	  cerr << "A" << "\n";
	How_Many_Things_To_Send_I(World,counts_out,counts_in);
	Send_Data_Somewhere_No_Block(World,counts_out,counts_in,integers,doubles,
				     dataI_out,dataI_in_send,how_manyI,
				     dataR_out,dataR_in_send,how_manyR);
	//	if(Rank == 0)
	//	  cerr << "AA" << "\n";
      }
    else if(World == HypreWorld)
      {
	//	if(Rank == 0)
	//	  cerr << "B" << "\n";
	Send_Data_Hypre_Directions(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in_send,how_manyI,
				   dataR_out,dataR_in_send,how_manyR);
	//	if(Rank == 0)
	//	  cerr << "BB" << "\n";
      }
    else if(tag==0 || tag==4)
      {
	//	if(Rank == 0)
	//	  cerr << "C" << "\n";
	Send_Data_Other_Directions(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in_send,how_manyI,
				   dataR_out,dataR_in_send,how_manyR);
	//	if(Rank == 0)
	//	  cerr << "CC" << "\n";
      }
    else
      {
	//	if(Rank == 0)
	//	  cerr << "D" << "\n";
	Send_Data_One_Directions(counts_out,counts_in,integers,doubles,
				 dataI_out,dataI_in_send,how_manyI,
				 dataR_out,dataR_in_send,how_manyR);
	//	if(Rank == 0)
	//	  cerr << "DD" << "\n";
      }
  }
  //
  void Mess::Send_Data_One_Directions(vector <int>& counts_out,vector <int>& counts_in,int integers,int doubles,
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
    //      FF << " SendOne AA " << FractalRank << "\n";
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
	cerr << " DUMP IT A " << FractalRank << " " << ba.what() << " " << totals << " " << nIdata << " " << nRdata << "\n";
	FF << " DUMP IT A " << ba.what() << " " << totals << " " << nIdata << " " << nRdata << "\n";
	for(int FR2=0;FR2<FractalNodes2;FR2++)
	  FF << FR2 << " " << dataI_out[FR2].size() << " " << dataR_out[FR2].size() << "\n";
	FF << endl;
	assert(0);
      }
    //      FF << " SendOne BB " << FractalRank << "\n";
    How_Many_Things_To_Send_I(MComms[2],countsa_out,countsa_in);
    int total_in=0;
    int total_out=0;
    for(int FR2=0;FR2<FractalNodes2;FR2++)
      {
	total_out+=countsa_out[FR2];
	total_in+=countsa_in[FR2];
      }
    //      FF << " TotalsOne 2 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << "\n";
    Send_Data_Somewhere_No_Block(MComms[2],countsa_out,countsa_in,
				 integers+1,doubles,
				 dataI_out,dataI_in,how_manyI,
				 dataR_out,dataR_in,how_manyR);
    //      FF << " SendOne CC " << FractalRank << "\n";
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
	cerr << " DUMP IT C " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT C " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes2;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    Full_Stop_Do_Not_Argue(MComms[2]);
    Full_Stop_Do_Not_Argue(MComms[1]);
    //      FF << " SendOne DD " << FractalRank << "\n";
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
    //      FF << " TotalsOne 1 " << total_out*(integers+2) << " " << total_out*doubles << " " <<  total_in*(integers+2) << " " << total_in*doubles << "\n";
    Send_Data_Somewhere_No_Block(MComms[1],countsa_out,countsa_in,
				 integers+2,doubles,
				 dataI_out,dataI_in,how_manyI,
				 dataR_out,dataR_in,how_manyR);

    //      FF << " SendOne EE " << FractalRank << "\n";
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
	cerr << " DUMP IT D " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT D " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes1;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    Full_Stop_Do_Not_Argue(MComms[1]);
    Full_Stop_Do_Not_Argue(MComms[0]);
    //      FF << " SendOne FF " << FractalRank << "\n";
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
    //      FF << " TotalsOne 0 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << "\n";
    Send_Data_Somewhere_No_Block(MComms[0],countsa_out,countsa_in,
				 integers+1,doubles,
				 dataI_out,dataI_in,how_manyI,
				 dataR_out,dataR_in,how_manyR);

    //      FF << " SendOne GG " << FractalRank << "\n";
    dataI_out.clear();
    dataR_out.clear();
    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    counts_in.assign(FractalNodes,0);
    countI=0;
    countR=0;
    try
      {
	//	  cerr << " testtestA " << FractalNodes0 << " " << integers << " " << doubles << "\n";
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
	//	  cerr << "testtestB " << FractalNodes0 << " " << integers << " " << doubles << endl;
      }
    catch(bad_alloc& ba)
      {
	cerr << " DUMP IT E " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT E " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes0;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    //      FF << "testtestC " << FractalNodes0 <<  " " << integers << " " << doubles << "\n";
    dataI_in.clear();
    dataR_in.clear();
    //      FF << " SendOne HH " << FractalRank <<  "\n";
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
	cerr << " DUMP IT F " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT F " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    Full_Stop_Do_Not_Argue();
    FF << " SendOne II " << FractalRank << " " << how_manyI << " " << how_manyR << "\n";
  }
  //
  void Mess::Send_Data_Other_Directions(vector <int>& counts_out,vector <int>& counts_in,int integers,int doubles,
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
    //      FF << " SendOther AA " << FractalRank << "\n";
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
	cerr << " DUMP IT A " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
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
	cerr << " DUMP IT B " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT B " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes0;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    dataI_in.clear();
    dataR_in.clear();
    //      FF << " SendOther BB " << FractalRank << "\n";
    How_Many_Things_To_Send_I(MComms[0],countsa_out,countsa_in);

    int total_in=0;
    int total_out=0;
    for(int FR0=0;FR0<FractalNodes0;FR0++)
      {
	total_out+=countsa_out[FR0];
	total_in+=countsa_in[FR0];
      }
    //      FF << " TotalsOther 0 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << "\n";

    Send_Data_Somewhere_No_Block(MComms[0],countsa_out,countsa_in,
				 integers+1,doubles,
				 dataI_out,dataI_in,how_manyI,
				 dataR_out,dataR_in,how_manyR);
    //      FF << " SendOther CC " << FractalRank << "\n";
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
	cerr << " DUMP IT C " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT C " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes1;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    Full_Stop_Do_Not_Argue(MComms[0]);
    Full_Stop_Do_Not_Argue(MComms[1]);
    //      FF << " SendOther DD " << FractalRank << "\n";
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
    //      FF << " TotalsOther 1 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << "\n";
    Send_Data_Somewhere_No_Block(MComms[1],countsa_out,countsa_in,
				 integers+2,doubles,
				 dataI_out,dataI_in,how_manyI,
				 dataR_out,dataR_in,how_manyR);

    //      FF << " SendOther EE " << FractalRank << "\n";
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
	cerr << " DUMP IT D " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT D " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes2;FR++)
	  FF << FR << " " << dataI_out[FR].size() << " " << dataR_out[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    Full_Stop_Do_Not_Argue(MComms[1]);
    Full_Stop_Do_Not_Argue(MComms[2]);
    //      FF << " SendOther FF " << FractalRank << "\n";
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
    //      FF << " TotalsOther 2 " << total_out*(integers+1) << " " << total_out*doubles << " " <<  total_in*(integers+1) << " " << total_in*doubles << "\n";
    Send_Data_Somewhere_No_Block(MComms[2],countsa_out,countsa_in,
				 integers+1,doubles,
				 dataI_out,dataI_in,how_manyI,
				 dataR_out,dataR_in,how_manyR);

    //      FF << " SendOther GG " << FractalRank << "\n";
    dataI_out.clear();
    dataR_out.clear();
    counts_in.assign(FractalNodes,0);
    countI=0;
    countR=0;
    counterI=0;
    counterR=0;
    vector < vector <int> > dataI;
    vector < vector <double> > dataR;
    try
      {
	for(int FR2=0;FR2<FractalNodes2;FR2++)
	  {
	    dataI.clear();
	    dataR.clear();
	    dataI.resize(FractalNodes01);
	    dataR.resize(FractalNodes01);
	    for(int c=0;c<countsa_in[FR2];c++)
	      {
		int FRFrom=dataI_in[countI];
		countI++;
		counts_in[FRFrom]++;
		int FR01=FRFrom % FractalNodes01;
		for(int nI=0;nI<integers;nI++)
		  {
		    dataI[FR01].push_back(dataI_in[countI]);
		    countI++;
		  }
		for(int nR=0;nR<doubles;nR++)
		  {
		    dataR[FR01].push_back(dataR_in[countR]);
		    countR++;
		  }
	      }
	    for(int FR01=0;FR01<FractalNodes01;FR01++)
	      {
		vector <int>::iterator pIb=dataI[FR01].begin();
		vector <int>::iterator pIe=dataI[FR01].end();
		while(pIb != pIe)
		  {
		    dataI_in[counterI]=*pIb;
		    pIb++;
		    counterI++;
		  }
		vector <double>::iterator pRb=dataR[FR01].begin();
		vector <double>::iterator pRe=dataR[FR01].end();
		while(pRb != pRe)
		  {
		    dataR_in[counterR]=*pRb;
		    pRb++;
		    counterR++;
		  }
	      }
	  }
	dataI_in.resize(counterI);
	dataR_in.resize(counterR);
	how_manyI=counterI;
	how_manyR=counterR;
      }
    catch(bad_alloc& ba)
      {
	cerr << " DUMP IT E " << FractalRank << " " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	FF << " DUMP IT E " << ba.what() << " " << dataI_in.size() << " " << dataR_in.size() << "\n";
	for(int FR=0;FR<FractalNodes01;FR++)
	  FF << FR << " " << dataI[FR].size() << " " << dataR[FR].size() << "\n";
	FF << endl;
	assert(0);
      }
    Full_Stop_Do_Not_Argue();
    FF << " SendOther II " << FractalRank << " " << how_manyI << " " << how_manyR << "\n";
  }
  //
  void Mess::Send_Data_Hypre_Directions(vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
					vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
					vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR)
  {
    double aNodes=HypreNodes;
    int HypreNodes0=pow(aNodes-0.5,1.0/3.0)+1.0;
    double a12=HypreNodes/HypreNodes0;
    int HypreNodes1=sqrt(a12-0.5)+1.0;
    int HypreNodes01=HypreNodes0*HypreNodes1;
    int HypreNodes2=HypreNodes/HypreNodes01;
    //    int HypreRank0=HypreRank % HypreNodes0;
    //    int HypreRank1=(HypreRank/HypreNodes0) % HypreNodes1;
    int HypreRank2=HypreRank/HypreNodes01;
    //      int HypreNodesBox=HypreNodes01*HypreNodes2;
    //      int extras=HypreNodes-HypreNodesBox;
    //      int ExtraLines=extras/HypreNodes0;
    //      int ExtraNodes=extras % HypreNodes0;
    int HypreLong2=how_many_nodes(HComms[2]);
    int HypreLong1=how_many_nodes(HComms[1]);
    int HypreLong0=how_many_nodes(HComms[0]);
//     if(HypreRank == 0)
//       {
// 	cerr << " AHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;
// 	cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
//       }
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
    //      cerr << " BHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    dataI_out.clear();
    dataR_out.clear();
    How_Many_Things_To_Send_I(HComms[2],countsa_out,countsa_in);
    //      cerr << " CHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";

    dataIa_in.clear();
    dataRa_in.clear();

    Send_Data_Somewhere_No_Block(HComms[2],countsa_out,countsa_in,
				 integers+1,doubles,
				 dataIa_out,dataIa_in,how_manyI,
				 dataRa_out,dataRa_in,how_manyR);
    //      cerr << " DHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
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
    //      cerr << " EHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    countsa_in.assign(HypreLong1,0);
    How_Many_Things_To_Send_I(HComms[1],countsa_out,countsa_in);
    dataIa_in.clear();
    dataRa_in.clear();
    //      cerr << " FHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";

    Send_Data_Somewhere_No_Block(HComms[1],countsa_out,countsa_in,
				 integers+2,doubles,
				 dataIa_out,dataIa_in,how_manyI,
				 dataRa_out,dataRa_in,how_manyR);
    //      cerr << " GHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";

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
    //      cerr << " HHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    Full_Stop_Do_Not_Argue(HComms[1]);
    Full_Stop_Do_Not_Argue(HComms[0]);
    countsa_in.assign(HypreLong0,0);
    How_Many_Things_To_Send_I(HComms[0],countsa_out,countsa_in);
    //      cerr << " IHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    dataIa_in.clear();
    dataRa_in.clear();
    Send_Data_Somewhere_No_Block(HComms[0],countsa_out,countsa_in,
				 integers+2,doubles,
				 dataIa_out,dataIa_in,how_manyI,
				 dataRa_out,dataRa_in,how_manyR);
    //      cerr << " JHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
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
    //      cerr << " KHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    Full_Stop_Do_Not_Argue(HComms[0]);
    Full_Stop_Do_Not_Argue(HComms[2]);
    countsa_in.assign(HypreLong2,0);
    How_Many_Things_To_Send_I(HComms[2],countsa_out,countsa_in);
    //      cerr << " LHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
    dataIa_in.clear();
    dataRa_in.clear();
    Send_Data_Somewhere_No_Block(HComms[2],countsa_out,countsa_in,
				 integers+1,doubles,
				 dataIa_out,dataIa_in,how_manyI,
				 dataRa_out,dataRa_in,how_manyR);
    //      cerr << " MHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
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
    //      cerr << " KHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
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
    //      cerr << " LHyp " << FractalRank << " " << HypreRank << " "  << HypreRank0 << " "  << HypreRank1 << " "  << HypreRank2;;
    //      cerr << " " << HypreLong0 << " " << HypreLong1 << " " << HypreLong2 << " " << HypreNodes << "\n";
  }
  void Mess::MPI_MYTest(int which,int test) const
  {
    if(test == MPI_SUCCESS)
      return;
    fprintf(p_file->PFFractalMemory," MPI Error %d %d %d %d %d %d %d %d \n",which,test,
	    MPI_ERR_COMM,MPI_ERR_TYPE,MPI_ERR_COUNT,MPI_ERR_TAG,MPI_ERR_RANK,MPI_SUCCESS);
  }
  void Mess::my_AllgatherI(vector <int>& paramsend,vector <int>& paramrecv,const int& nsend)
  {
    int ROOT=ROOTNODE;
    MPI_Gather(&(*(paramsend.begin())),nsend,MPI_INT,&(*(paramrecv.begin())),nsend,MPI_INT,ROOT,FractalWorld);
    MPI_Bcast(&(*paramrecv.begin()),nsend*FractalNodes,MPI_INT,ROOT,FractalWorld);
  }
  void Mess::my_AllgatherI(vector <long int>& paramsend,vector <long int>& paramrecv,const int& nsend)
  {
    int ROOT=ROOTNODE;
    MPI_Gather(&(*(paramsend.begin())),nsend,MPI_LONG,&(*(paramrecv.begin())),nsend,MPI_LONG,ROOT,FractalWorld);
    MPI_Bcast(&(*paramrecv.begin()),nsend*FractalNodes,MPI_LONG,ROOT,FractalWorld);
  }
  void Mess::my_AllgatherR(vector <double>& paramsend,vector <double>& paramrecv,const int& nsend) const
  {
    int ROOT=ROOTNODE;
    MPI_Gather(&(*(paramsend.begin())),nsend,MPI_DOUBLE,&(*(paramrecv.begin())),nsend,MPI_DOUBLE,ROOT,FractalWorld);
    Send_DOUBLE_from_ROOT(paramrecv,nsend*FractalNodes,ROOT);
  }
  void Mess::calc_total_particles(const int& NP)
  {
    vector <long int> particles(1);
    particles[0]=NP;
    Find_Sum_LONG_INT(particles,1);
    number_particles_total=particles[0];
  }
  void Mess::Find_Max_INT(vector <int>& integers,const int& how_long) const
  {
    int ROOT=ROOTNODE;
    Find_Max_INT_to_ROOT(integers,how_long,ROOT);
    Send_INT_from_ROOT(integers,how_long,ROOT);
  }
  void Mess::Find_Max_DOUBLE(vector <double>& doubles,const int& how_long) const
  {
    int ROOT=ROOTNODE;
    Find_Max_DOUBLE_to_ROOT(doubles,how_long,ROOT);
    Send_DOUBLE_from_ROOT(doubles,how_long,ROOT);
  }
  void Mess::Find_Sum_LONG_INT(vector <long int>& integers,const int& how_long) const
  {
    int ROOT=ROOTNODE;
    Find_Sum_LONG_INT_to_ROOT(integers,how_long,ROOT);
    Send_LONG_INT_from_ROOT(integers,how_long,ROOT);
  }
  void Mess::Find_Sum_INT(vector <int>& integers,const int& how_long) const
  {
    int ROOT=ROOTNODE;
    Find_Sum_INT_to_ROOT(integers,how_long,ROOT);
    Send_INT_from_ROOT(integers,how_long,ROOT);
  }
  void Mess::Find_Sum_DOUBLE(vector <double>& doubles,const int& how_long) const
  {
    int ROOT=ROOTNODE;
    Find_Sum_DOUBLE_to_ROOT(doubles,how_long,ROOT);
    Send_DOUBLE_from_ROOT(doubles,how_long,ROOT);
  }
  void Mess::Send_INT_from_ROOT(int* numbers,const int& how_long,const int& ROOT) const
  {
    MPI_Bcast(numbers,how_long,MPI_INT,ROOT,FractalWorld);
  }
  void Mess::Find_Max_INT_to_ROOT(vector <int>& numbers,const int& how_long,const int& ROOT) const
  {
    vector <int> maxi(how_long);
    MPI_Reduce(&(*numbers.begin()),&(*maxi.begin()),how_long,MPI_INT,MPI_MAX,ROOT,FractalWorld);
    numbers=maxi;
  }
  void Mess::Find_Max_DOUBLE_to_ROOT(vector <double>& numbers,const int& how_long,const int& ROOT) const
  {
    vector <double> maxr(how_long);
    MPI_Reduce(&(*numbers.begin()),&(*maxr.begin()),how_long,MPI_DOUBLE,MPI_MAX,ROOT,FractalWorld);
    numbers=maxr;
  }
  void Mess::Find_Sum_INT_to_ROOT(vector <int>& numbers,const int& how_long,const int& ROOT) const
  {
    vector <int> sumup(how_long);
    MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_INT,MPI_SUM,ROOT,FractalWorld);
    numbers=sumup;
  }
  void Mess::Find_Sum_LONG_INT_to_ROOT(vector <long int>& numbers,const int& how_long,const int& ROOT) const
  {
    vector <long int> sumup(how_long);
    MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_LONG,MPI_SUM,ROOT,FractalWorld);
    numbers=sumup;
  }
  //   void Mess::Find_Sum_DOUBLE_to_ROOT(vector < vector < vector <double> > >& numbers,const int& how_long,vector <int>& ROOTS) const
  //   {
  //     vector <double> sumup(how_long);
  //     vector <MPI_Request>Mreq;
  //     int FR=0;
  //     for(int FR2=0;FR2<FractalNodes2;FR2++)
  //       for(int FR1=0;FR1<FractalNodes1;FR1++)
  // 	{
  // 	  Mreq.push_back(MPI_Request());
  // 	  MPI_Ireduce(&(*numbers[FR2][FR1].begin()),&(*sumup.begin()),how_long,MPI_DOUBLE,MPI_SUM,ROOTS[FR],FractalWorld,Mreq.back());
  // 	}
  //     vector <MPI_Status> Mstat(FractalNodes2*FractalNodes1);
  //     int answer=MPI_Waitall(FractalNodes2*FractalNodes1,&(*Mreq.begin()),&(*Mstat.begin()));
  //     MPI_MYTest(34,answer);
  //     FR=0;
  //     for(int FR2=0;FR2<FractalNodes2;FR2++)
  //       for(int FR1=0;FR1<FractalNodes1;FR1++)
  // 	{
  // 	  if(ROOTS[FR] == FractalRank)
  // 	    numbers[FR2][FR1]=sumup;
  // 	  FR++;
  // 	}
  //   }
  //   void Mess::Find_Sum_DOUBLE_to_ROOT(vector < vector <double> >& numbers,const int& how_long,vector <int>& ROOTS) const
  //   {
  //     vector <double> sumup(how_long);
  //     vector <MPI_Request>Mreq;
  //     for(int FR2=0;FR2<FractalNodes2;FR2++)
  //       {
  // 	Mreq.push_back(MPI_Request());
  // 	MPI_Ireduce(&(*numbers[FR2].begin()),&(*sumup.begin()),how_long,MPI_DOUBLE,MPI_SUM,ROOTS[FR2],FractalWorld,Mreq.back());
  //       }
  //     vector <MPI_Status> Mstat(FractalNodes2);
  //     int answer=MPI_Waitall(FractalNodes2,&(*Mreq.begin()),&(*Mstat.begin()));
  //     MPI_MYTest(14,answer);
  //     for(int FR2=0;FR2<FractalNodes2;FR2++)
  //       if(ROOTS[FR2] == FractalRank)
  // 	numbers[FR2]=sumup;
  //   }
  void Mess::Find_Sum_DOUBLE_to_ROOT(vector <double>& numbers,const int& how_long,const int& ROOT) const
  {
    vector <double> sumup(how_long);
    MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_DOUBLE,MPI_SUM,ROOT,FractalWorld);
    numbers=sumup;
  }
  void Mess::Find_Sum_FLOAT_to_ROOT(vector <float>& numbers,const int& how_long,const int& ROOT,MPI_Comm& World) const
  {
    vector <float> sumup(how_long);
    MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_FLOAT,MPI_SUM,ROOT,World);
    numbers=sumup;
  }
  void Mess::Find_Sum_FLOAT_to_ROOT(vector <float>& numbers,const int& how_long,const int& ROOT) const
  {
    vector <float> sumup(how_long);
    if(FractalNodes  <= MPI_SWITCH)
      {
	MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_FLOAT,MPI_SUM,ROOT,FractalWorld);
	numbers=sumup;
	return;
      }
    int ROOT2=ROOT/(FractalNodes0*FractalNodes1);
    int ROOT1=(ROOT/FractalNodes0) % FractalNodes1;
    int ROOT0=ROOT % FractalNodes0;
    //      cerr << " Reduce " << FractalRank << " " << ROOT << " " << ROOT0 << " " << ROOT1 << " " << ROOT2 ;
    //      cerr << " " << how_long << " " << numbers.size() << " " << sumup.size() << " " << MComms.size() << endl;
    MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_FLOAT,MPI_SUM,ROOT2,MComms[2]);
    numbers=sumup;
    Full_Stop_Do_Not_Argue();
    if(ROOT2 == FractalRank2)
      {
	MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_FLOAT,MPI_SUM,ROOT1,MComms[1]);
	numbers=sumup;
      }
    Full_Stop_Do_Not_Argue();
    if(ROOT1 == FractalRank1 && ROOT2 == FractalRank2)
      {
	MPI_Reduce(&(*numbers.begin()),&(*sumup.begin()),how_long,MPI_FLOAT,MPI_SUM,ROOT0,MComms[0]);
	numbers=sumup;
      }
    Full_Stop_Do_Not_Argue();
  }
  //   void Mess::Send_INT_from_ROOT(vector < vector < vector <int> > >& numbers,const int& how_long,const vector <int>& ROOTS) const
  //   {
  //     vector <MPI_Request> Mreq;
  //     int FR=0;
  //     for(int FR2=0;FR2<FractalNodes2;FR2++)
  //       for(int FR1=0;FR1<FractalNodes1;FR1++)
  // 	{
  // 	  Mreq.push_back(MPI_Request());
  // 	  MPI_Ibcast(&(*numbers[FR2][FR1].begin()),how_long,MPI_INT,ROOTS[FR],FractalWorld,Mreq.back());
  // 	  FR++;
  // 	}
  //     vector <MPI_Status> Mstat(FractalNodes2*FractalNodes1);
  //     int answer=MPI_Waitall(FractalNodes2*FractalNodes1,&(*Mreq.begin()),&(*Mstat.begin()));
  //     MPI_MYTest(24,answer);
  //   }
  //   void Mess::Send_INT_from_ROOT(vector < vector <int> >& numbers,const int& how_long,const vector <int>& ROOTS) const
  //   {
  //     vector <MPI_Request> Mreq;
  //     for(int FR2=0;FR2<FractalNodes2;FR2++)
  //       {
  // 	Mreq.push_back(MPI_Request());
  // 	MPI_Ibcast(&(*numbers[FR2].begin()),how_long,MPI_INT,ROOTS[FR2],FractalWorld,Mreq.back());
  //       }
  //     vector <MPI_Status> Mstat(FractalNodes2);
  //     int answer=MPI_Waitall(FractalNodes2,&(*Mreq.begin()),&(*Mstat.begin()));
  //     MPI_MYTest(14,answer);
  //   }
  void Mess::Send_INT_from_ROOT(vector <int>& numbers,const int& how_long,const int& ROOT) const
  {
    MPI_Bcast(&(*numbers.begin()),how_long,MPI_INT,ROOT,FractalWorld);
  }
  void Mess::Send_LONG_INT_from_ROOT(vector <long int>& numbers,const int& how_long,const int& ROOT) const
  {
    MPI_Bcast(&(*numbers.begin()),how_long,MPI_LONG,ROOT,FractalWorld);
  }
  void Mess::Send_DOUBLE_from_ROOT(vector <double>& numbers,const int& how_long,const int& ROOT) const
  {
    MPI_Bcast(&(*numbers.begin()),how_long,MPI_DOUBLE,ROOT,FractalWorld);
  }
  void Mess::Full_Stop() const
  {
    if(time_trial)
      MPI_Barrier(FractalWorld);
  }
  void Mess::Full_Stop_Do_Not_Argue() const
  {
    MPI_Barrier(FractalWorld);
  }
  void Mess::Full_Stop(MPI_Comm& World) const
  {
    if(time_trial)
      MPI_Barrier(World);
  }
  void Mess::Full_Stop_Do_Not_Argue(MPI_Comm& World) const
  {
    MPI_Barrier(World);
  }
  void Mess::zeroR()
  {
    zeroR(0.0);
  }
  void Mess::zeroR(const double& grail)
  {
    std::fill(potR,potR+2*total_memory,grail);
  }
  void Mess::zeroRS()
  {
    zeroRS(0.0);
  }
  void Mess::zeroRS(const double& grail)
  {
    if(IAmPeriodic)
      std::fill(potRS,potRS+2*total_memory,grail);
    else
      std::fill(potRS,potRS+(glength+1)*(glength+1)*length_x,grail);
  }
  double Mess::Clock() const
  {
    return MPI_Wtime();
  }
  int Mess::MyHypreRank() const
  {
    int rank;
    if(HypreWorld == MPI_COMM_NULL)
      return -1;
    MPI_Comm_rank(HypreWorld,&rank);
    return rank;      
  }
  void Mess::HypreGroupCreate(vector <int>& ranks)
  {
    MPI_Comm_group(FractalWorld,&FractalGroup);
    MPI_Group_incl(FractalGroup, HypreNodes, &(*(ranks.begin())), &HypreGroup);
    MPI_Comm_create(FractalWorld, HypreGroup, &HypreWorld);
    HypreRank=MyHypreRank();
    if(!IAmAHypreNode)
      return;
    make_MPI_Hypre_Groups();
  }
  void Mess::HypreGroupFree()
  {
    if(!IAmAHypreNode)
      return;
    MPI_Group_free(&HypreGroup);
    MPI_Comm_free(&HypreWorld);
    HypreWorld=MPI_COMM_NULL;
    HypreGroups3Free();
  }
  void Mess::HypreGroups3Free()
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
  void Mess::createFractalWorld(MPI_Comm& World,vector <int>& dims)
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
  void Mess::make_Random_Group()
  {
    if(FractalNodes <= 2*RandomNodes)
      {
	RandomNodes=FractalNodes;
	RandomWorld=FractalWorld;
	RandomGroup=FractalGroup;
	IAmARandomNode=true;
	Rranks.resize(FractalNodes);
	IRranks.resize(FractalNodes);
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    Rranks[FR]=FR;
	    IRranks[FR]=FR;
	  }
      }
    else
      {
	IAmARandomNode=false;
	Rranks.assign(RandomNodes,-1);
	IRranks.assign(FractalNodes,-1);
	double aFractalNodes=FractalNodes;
	double rand_max=RAND_MAX;
	for(int RR=0;RR<RandomNodes;RR++)
	  {
	    bool not_yet=true;
	    while(not_yet)
	      {
		int FR=aFractalNodes*(double)(rand())/rand_max;
		if(IRranks[FR] < 0)
		  {
		    Rranks[RR]=FR;
		    IRranks[FR]=RR;
		    not_yet=false;
		  }
	      }
	  }
	std::sort(Rranks.begin(),Rranks.end());
	IRranks.assign(FractalNodes,-1);
	for(int RR=0;RR<RandomNodes;RR++)
	  {
	    IRranks[Rranks[RR]]=RR;
	    if(Rranks[RR] == FractalRank)
	      {
		RandomRank=RR;
		IAmARandomNode=true;
	      }
	  }
      }
    MPI_Comm_group(FractalWorld,&FractalGroup);
    MPI_Group_incl(FractalGroup, RandomNodes, &(*(Rranks.begin())), &RandomGroup);
    MPI_Comm_create(FractalWorld, RandomGroup, &RandomWorld);
    if(!IAmARandomNode)
      return;
    int Ranky;
    MPI_Comm_rank(RandomWorld,&Ranky);
    assert(Ranky == RandomRank);
  }
}
