#ifndef _Mess_Defined_
#define _Mess_Defined_

#include "fftw3-mpi.h"

#include <fstream>
#include <vector>
#include <deque>

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
    int WallNodes;
    const int N63;
    int ROOTNODE;
    int FFTRank;
    int FFTNodes;
    int HypreRank;
    int HypreNodes;
    int mynumber;
    int MPI_SWITCH;
    int MPI_MAX_COMMS;
    long int number_particles_total;
    std::deque <Particle*> parts_tmp;
    std::deque <Particle*> parts_tmpp;
    Particle* Parts_in;
    std::vector < std::vector <int> > Slices;
    std::vector < std::vector <int> > BoxS;
    std::vector < std::vector <int> > BoxSL;
    std::vector < std::vector <bool> > counts_on_nodes;
    std::vector <bool> count_on_node;
    std::vector<std::vector <int>> node_lists;
    std::vector <int> freenodes;
    int glength;
    pint start_x;
    pint length_x;
    pint total_memory;
    fftw_plan plan_rc;
    fftw_plan plan_cr;
    double* potR;
    double* potRS;
    fftw_complex* potC; 
    std::vector <int> WhichSlice;
    std::vector <double> green;
    MPI_Comm FractalWorld;
    MPI_Comm FFTWorld;
    MPI_Comm HypreWorld;
    MPI_Group FractalGroup;
    MPI_Group FFTGroup;
    MPI_Group HypreGroup;
    bool IAmPeriodic;
    bool IAmAnFFTNode;
    bool IAmAHypreNode;
    std::vector <int>Franks;
    std::vector <bool>ItIsAnFFTNode;
    std::vector <int>Hranks;
    std::vector <int>IHranks;
    std::vector <int>Rranks;
    std::vector <int>IRranks;
    std::vector <MPI_Comm> MComms;
    std::vector <MPI_Comm> HComms;
    std::vector <MPI_Group> HG;
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
      N63(-1),
      ROOTNODE(0),
      FFTRank(0),
      FFTNodes(1234567),
      HypreRank(0),
      HypreNodes(0),
      mynumber(-1),
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
      std::cerr << " Empty Mess " << "\n";
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP,
	 int& FR0,int& FR1,int& FR2,const int& FN,MPI_Comm& FW):
      FractalRank(0),
      FractalNodes(1),
      FractalNodes0(FR0),
      FractalNodes1(FR1),
      FractalNodes2(FR2),
      N63(-1),
      ROOTNODE(0),
      FFTNodes(FN),
      HypreRank(0),
      HypreNodes(0),
      mynumber(-1),
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
      int grid_length=GR;
      IAmPeriodic=PR;
      bool periodic=PR;
      WallTime=Clock();
      if(MR)
	{
	  WallNodes=FractalNodes0*FractalNodes1*FractalNodes2-
	    (FractalNodes0-2)*(FractalNodes1-2)*(FractalNodes2-2);
	  MPIStartup(PR,FR0,FR1,FR2);
	  FractalRank=what_is_my_rank(); 
	  FractalNodes=how_many_nodes(); 
	  assert(FractalNodes == FR0*FR1*FR2);
	  FractalRank0=FractalRank % FractalNodes0;
	  FractalRank1=(FractalRank/FractalNodes0) % FractalNodes1;
	  FractalRank2=FractalRank/(FractalNodes0*FractalNodes1);
	  FFTWStartup(grid_length,periodic);
	  ROOTNODE=0;
	  calc_fftw_Slices(grid_length,periodic);	
	  calc_total_particles(NP);
	  make_MPI_Groups();
	}
      else
	{
	  number_particles_total=NP;
	  length_x=grid_length;
	}
    }
    Mess(const bool& MR,const int& GR,const bool& PR,const int& NP,const int& FN,MPI_Comm& FW):
      FractalRank(0),
      FractalNodes(1),
      ROOTNODE(0),
      N63(-1),
      FFTNodes(FN),
      HypreRank(0),
      HypreNodes(0),
      mynumber(-1),
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
    }
    ~Mess()
    {
      FFTWFinal();
      if(standalone)
	MPIFinal();
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
    void dumpR(std::ofstream& FILE,const int& length) const;
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
    template <class T> void How_Many_On_Nodes(T count,std::vector <T>& counts) const;
    void MAX_Things_To_Send_Receive_I(std::vector <int>& counts_out_send,std::vector <int>& counts_in_send,std::vector <int>& maxSR);
    long int How_Many_In_Solver(const int S) const;
    void How_Many_Things_To_Send_I(std::vector <int>& counts_out_send,std::vector <int>& counts_in_send);
    void How_Many_Things_To_Send_I(MPI_Comm& World,
					 std::vector <int>& counts_out_send,std::vector <int>& counts_in_send);
    void Send_Data_Somewhere_No_Block(std::vector <int>& counts_out_send,std::vector <int>& counts_in_send,const int& integers,const int& doubles,
					    std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in_send,int& how_manyI,
					    std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in_send,int& how_manyR);
    void Send_Data_Somewhere_No_Block(MPI_Comm& World,
					    std::vector <int>& counts_out_send,std::vector <int>& counts_in_send,const int& integers,const int& doubles,
					    std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in_send,int& how_manyI,
					    std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in_send,int& how_manyR);
    void make_MPI_Hypre_Groups();
    void make_MPI_Groups();
    void Send_Data_Some_How(int tag,std::vector <int>& counts_out,std::vector <int>& counts_in,const int& integers,const int& doubles,
				  std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in_send,int& how_manyI,
				  std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in_send,int& how_manyR);
    void Send_Data_Some_How(int tag,MPI_Comm& World,
				  std::vector <int>& counts_out,std::vector <int>& counts_in,const int& integers,const int& doubles,
				  std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in_send,int& how_manyI,
				  std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in_send,int& how_manyR);
    void Send_Data_One_Directions(std::vector <int>& counts_out,std::vector <int>& counts_in,int integers,int doubles,
					std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in,int& how_manyI,
					std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in,int& how_manyR);
    void Send_Data_Other_Directions(std::vector <int>& counts_out,std::vector <int>& counts_in,int integers,int doubles,
					  std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in,int& how_manyI,
					  std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in,int& how_manyR);
    void Send_Data_Hypre_Directions(std::vector <int>& counts_out,std::vector <int>& counts_in,const int& integers,const int& doubles,
					  std::vector < std::vector <int> >& dataI_out,std::vector <int>& dataI_in,int& how_manyI,
					  std::vector < std::vector <double> >& dataR_out,std::vector <double>& dataR_in,int& how_manyR);
    template <class GO_AWAY> void really_clear(std::vector <GO_AWAY>& die);
    void MPI_MYTest(int which,int test) const;
    void Which_Nodes(int count,std::vector <int>& counts,std::vector <bool>& YesNo,int ROOT,MPI_Comm& World);
    void my_AllgatherI(MPI_Comm& World,std::vector <int>& paramsend,std::vector <int>& paramrecv,const int& nsend) const;
    void my_AllgatherI(std::vector <int>& paramsend,std::vector <int>& paramrecv,const int& nsend) const;
    void my_AllgatherI(std::vector <long>& paramsend,std::vector <long>& paramrecv,const int& nsend) const;
    void my_AllgatherR(std::vector <double>& paramsend,std::vector <double>& paramrecv,const int& nsend) const;
    void calc_total_particles(const int& NP);
    void Find_Max_INT(MPI_Comm& World,std::vector <int>& integers,const int& how_long) const;
    void Find_Max_INT(std::vector <int>& integers,const int& how_long) const;
    void Find_Max_DOUBLE(std::vector <double>& doubles,const int& how_long) const;
    void Find_Sum_LONG_INT(std::vector <long int>& integers,const int& how_long) const;
    void Find_Sum_INT(std::vector <int>& integers,const int& how_long) const;
    void Find_Sum_DOUBLE(std::vector <double>& doubles,const int& how_long) const;
    void Send_INT_from_ROOT(MPI_Comm& World,int* numbers,const int& how_long,const int& ROOT) const;
    void Send_INT_from_ROOT(int* numbers,const int& how_long,const int& ROOT) const;
    void Find_Max_INT_to_ROOT(MPI_Comm& World,std::vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Max_INT_to_ROOT(std::vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Max_DOUBLE_to_ROOT(std::vector <double>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_FLOAT_to_ROOT(std::vector <float>& numbers,const int& how_long,const int& ROOT,MPI_Comm& World) const;
    void Find_Sum_FLOAT_to_ROOT(std::vector <float>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_INT_to_ROOT(std::vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_LONG_INT_to_ROOT(std::vector <long int>& numbers,const int& how_long,const int& ROOT) const;
    void Find_Sum_DOUBLE_to_ROOT(std::vector < std::vector < std::vector <double> > >& numbers,const int& how_long,std::vector <int>& ROOTS) const;
    void Find_Sum_DOUBLE_to_ROOT(std::vector < std::vector <double> >& numbers,const int& how_long,std::vector <int>& ROOTS) const;
    void Find_Sum_DOUBLE_to_ROOT(std::vector <double>& numbers,const int& how_long,const int& ROOT) const;
    void Send_INT_from_ROOT(std::vector < std::vector < std::vector <int> > >& numbers,const int& how_long,const std::vector <int>& ROOTS) const;
    void Send_INT_from_ROOT(std::vector < std::vector <int> >& numbers,const int& how_long,const std::vector <int>& ROOTS) const;
    void Send_INT_from_ROOT(std::vector <int>& numbers,const int& how_long,const int& ROOT) const;
    void Send_LONG_INT_from_ROOT(std::vector <long int>& numbers,const int& how_long,const int& ROOT) const;
    void Send_DOUBLE_from_ROOT(std::vector <double>& numbers,const int& how_long,const int& ROOT) const;
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
    void HypreGroupCreate(std::vector <int>& ranks);
    void HypreGroupFree();
    void HypreGroups3Free();
    void createFractalWorld(MPI_Comm& World,std::vector <int>& dims);
  };
}
#endif
