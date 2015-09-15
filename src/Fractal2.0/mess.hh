#ifndef _Mess_Defined_
#define _Mess_Defined_
namespace FractalSpace
{
  typedef ptrdiff_t pint;
  class Mess{
  public:
    int FractalRank;
    int FractalRank0;
    int FractalRank1;
    int FractalRank2;
    int FractalNodes;
    int FractalNodes0;
    int FractalNodes1;
    int FractalNodes2;
    int ROOTNODE;
    int FFTRank;
    int FFTNodes;
    int HypreRank;
    int HypreNodes;
//     int RandomRank;
//     int RandomNodes;
    int MPI_SWITCH;
    int MPI_MAX_COMMS;
    long int number_particles_total;
    vector <Particle*> parts_tmp;
    vector <Particle*> parts_tmpp;
    Particle* Parts_in;
    vector < vector <int> > Slices;
    vector < vector <int> > BoxS;
    vector < vector <int> > BoxSL;
    int glength;
    pint start_x;
    pint length_x;
    pint total_memory;
    fftw_plan plan_rc;
    fftw_plan plan_cr;
    double* potR;
    double* potRS;
    fftw_complex* potC; 
    vector <int> WhichSlice;
    vector <double> green;
    MPI_Comm FractalWorld;
    MPI_Comm FFTWorld;
    MPI_Comm HypreWorld;
//     MPI_Comm RandomWorld;
    MPI_Group FractalGroup;
    MPI_Group FFTGroup;
    MPI_Group HypreGroup;
//     MPI_Group RandomGroup;
    bool IAmPeriodic;
    bool IAmAnFFTNode;
    bool IAmAHypreNode;
//     bool IAmARandomNode;
    vector <int>Franks;
    vector <bool>ItIsAnFFTNode;
    vector <int>Hranks;
    vector <int>IHranks;
    vector <int>Rranks;
    vector <int>IRranks;
    vector <MPI_Comm> MComms;
    vector <MPI_Comm> HComms;
    vector <MPI_Group> HG;
    bool possibleDANGER;
    int DANGERlevel;
    int fftwTAG;
    bool time_trial;
    bool standalone;
    File* p_file;
    double WallTime;
    double TreeTime;
    static bool IAMROOT;
    Mess():
      FractalRank(0),
      FractalNodes(1),
      FFTRank(0),
      FFTNodes(1234567),
      HypreRank(0),
      HypreNodes(0),
      MPI_SWITCH(512),
      MPI_MAX_COMMS(4096),
      number_particles_total(-1),
      glength(1),
      start_x(0),
      length_x(-1),
      total_memory(-1),
      IAmPeriodic(true),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      standalone(true),
      TreeTime(-1.0)
    {
      WallTime=Clock();
      cerr << " Empty Mess " << "\n";
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
      MPI_MAX_COMMS(4096),
      number_particles_total(-1),
      glength(GR),
      start_x(0),
      length_x(GR),
      total_memory(1),
      FractalWorld(FW),
      IAmPeriodic(PR),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      TreeTime(-1.0)
    {
      //      cerr << " Making a Mess with parameters" << "\n";
      int grid_length=GR;
      IAmPeriodic=PR;
      bool periodic=PR;
      WallTime=Clock();
      if(MR)
	{
	  MPIStartup(PR,FR0,FR1,FR2);
	  FractalRank=what_is_my_rank(); 
	  FractalNodes=how_many_nodes(); 
	  assert(FractalNodes == FR0*FR1*FR2);
	  FractalRank0=FractalRank % FractalNodes0;
	  FractalRank1=(FractalRank/FractalNodes0) % FractalNodes1;
	  FractalRank2=FractalRank/(FractalNodes0*FractalNodes1);
	  ROOTNODE=(FractalNodes0+FractalNodes0*FractalNodes1+FractalNodes)/2;
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
      //      if(FractalRank == 0)
      //	cerr << " made a mess " << FractalRank << " " << FractalNodes << " " << length_x << " " << start_x << " " << total_memory << "\n";
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP,const int& FN,MPI_Comm& FW):
      FractalRank(0),
      FractalNodes(1),
      FFTNodes(FN),
      HypreRank(0),
      HypreNodes(0),
      MPI_SWITCH(512),
      MPI_MAX_COMMS(4096),
      number_particles_total(-1),
      glength(GR),
      start_x(0),
      length_x(GR),
      total_memory(1),
      FractalWorld(FW),
      IAmPeriodic(PR),
      IAmAnFFTNode(true),
      IAmAHypreNode(true),
      time_trial(true),
      TreeTime(-1.0)
    {
      //      cerr << " Making a Mess with parameters" << "\n";
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
      //      if(FractalRank == 0)
      //	cerr << " made a mess " << FractalRank << " " << FractalNodes << " " << length_x << " " << start_x << " " << total_memory << "\n";
    }
    ~Mess()
    {
      //      cerr << " starting to clean up a mess " << FractalRank << "\n";
      FFTWFinal();
      if(standalone)
	MPIFinal();
      //      cerr << " cleaned up a mess " << FractalRank << "\n";
    }
    void MPIStartup();
    void MPIStartup(const bool& PR,int& FR0,int& FR1,int& FR2);
    void MPIFinal() const;
    int what_is_my_rank() const;
    int what_is_my_rank(MPI_Comm& World) const;
    int how_many_nodes(MPI_Comm& World) const;
    int how_many_nodes() const;
    int what_is_my_FFT_rank() const;
    int how_many_FFT_nodes() const;
    int what_is_my_Hypre_rank() const;
    int how_many_Hypre_nodes() const;
    void FFTWStartup(const int& length_1,const bool& periodic);
    void doFFTWorld(int how_long,const bool& periodic);
    void FFTWFinal();
    void dumpR(ofstream& FILE,const int& length) const;
    void create_potRC();
    void free_potRC();
    void create_potR();
    void free_potR();
    void create_potRS();
    void free_potRS();
    void create_potC();
    void free_potC();
    void fftw_real_to_complex();
    void fftw_complex_to_real();
    int fftw_where(const int& i,const int& j,const int& k,const int& lb,const int& lc) const;
    void calc_fftw_Slices(const int& length_a,const bool& periodic);
    template <class T> void How_Many_On_Nodes(T count,vector <T>& counts) const;
    void MAX_Things_To_Send_Receive_I(vector <int>& counts_out_send,vector <int>& counts_in_send,vector <int>& maxSR);
    void How_Many_Things_To_Send_I(vector <int>& counts_out_send,vector <int>& counts_in_send);
    void How_Many_Things_To_Send_I(MPI_Comm& World,
					 vector <int>& counts_out_send,vector <int>& counts_in_send);
    void Send_Data_Somewhere_No_Block(vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
					    vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
					    vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR);
    void Send_Data_Somewhere_No_Block(MPI_Comm& World,
					    vector <int>& counts_out_send,vector <int>& counts_in_send,const int& integers,const int& doubles,
					    vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
					    vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR);
    void make_MPI_Hypre_Groups();
    void make_MPI_Groups();
    void Send_Data_Some_How(int tag,vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
				  vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
				  vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR);
    void Send_Data_Some_How(int tag,MPI_Comm& World,
				  vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
				  vector < vector <int> >& dataI_out,vector <int>& dataI_in_send,int& how_manyI,
				  vector < vector <double> >& dataR_out,vector <double>& dataR_in_send,int& how_manyR);
    void Send_Data_One_Directions(vector <int>& counts_out,vector <int>& counts_in,int integers,int doubles,
					vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
					vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR);
    void Send_Data_Other_Directions(vector <int>& counts_out,vector <int>& counts_in,int integers,int doubles,
					  vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
					  vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR);
    void Send_Data_Hypre_Directions(vector <int>& counts_out,vector <int>& counts_in,const int& integers,const int& doubles,
					  vector < vector <int> >& dataI_out,vector <int>& dataI_in,int& how_manyI,
					  vector < vector <double> >& dataR_out,vector <double>& dataR_in,int& how_manyR);
    void MPI_MYTest(int which,int test) const;
    void Which_Nodes(int count,vector <int>& counts,vector <bool>& YesNo,int ROOT,MPI_Comm& World);
    void my_AllgatherI(vector <int>& paramsend,vector <int>& paramrecv,const int& nsend) const;
    void my_AllgatherI(vector <long>& paramsend,vector <long>& paramrecv,const int& nsend) const;
    void my_AllgatherR(vector <double>& paramsend,vector <double>& paramrecv,const int& nsend) const;
    void calc_total_particles(const int& NP);
    void Find_Max_INT(vector <int>& integers,const int& how_long) const;
    void Find_Max_DOUBLE(vector <double>& doubles,const int& how_long) const;
    void Find_Sum_LONG_INT(vector <long int>& integers,const int& how_long) const;
    void Find_Sum_INT(vector <int>& integers,const int& how_long) const;
    void Find_Sum_DOUBLE(vector <double>& doubles,const int& how_long) const;
    void Send_INT_from_ROOT(int* numbers,const int& how_long,const int& ROOT) const;
    void Find_Max_INT_to_ROOT(vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Max_DOUBLE_to_ROOT(vector <double>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_FLOAT_to_ROOT(vector <float>& numbers,const int& how_long,const int& ROOT,MPI_Comm& World) const;
    void Find_Sum_FLOAT_to_ROOT(vector <float>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_INT_to_ROOT(vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_LONG_INT_to_ROOT(vector <long int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_DOUBLE_to_ROOT(vector < vector < vector <double> > >& numbers,const int& how_long,vector <int>& ROOTS) const;
    void Find_Sum_DOUBLE_to_ROOT(vector < vector <double> >& numbers,const int& how_long,vector <int>& ROOTS) const;
    void Find_Sum_DOUBLE_to_ROOT(vector <double>& numbers,const int& how_long,const int& ROOT) const;
    void Send_INT_from_ROOT(vector < vector < vector <int> > >& numbers,const int& how_long,const vector <int>& ROOTS) const;
    void Send_INT_from_ROOT(vector < vector <int> >& numbers,const int& how_long,const vector <int>& ROOTS) const;
    void Send_INT_from_ROOT(vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Send_LONG_INT_from_ROOT(vector <long int>& numbers,const int& how_long,const int& ROOT) const;
    void Send_DOUBLE_from_ROOT(vector <double>& numbers,const int& how_long,const int& ROOT) const;
    void Full_Stop() const;
    void Full_Stop(MPI_Comm& World) const;
    void Full_Stop_Do_Not_Argue() const;
    void Full_Stop_Do_Not_Argue(MPI_Comm& World) const;
    void Fake_Stop_Do_Not_Argue() const;
    void Fake_Stop_Do_Not_Argue(MPI_Comm& World) const;
    void zeroR();
    void zeroR(const double& grail);
    void zeroRS();
    void zeroRS(const double& grail);
    double Clock() const;
    int MyHypreRank() const;
    void HypreGroupCreate(vector <int>& ranks);
    void HypreGroupFree();
    void HypreGroups3Free();
    void createFractalWorld(MPI_Comm& World,vector <int>& dims);
//     void make_Random_Group();
  };
}
#endif
