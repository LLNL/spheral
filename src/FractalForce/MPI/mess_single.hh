#ifndef _Mess_Defined_
#define _Mess_Defined_
namespace FractalSpace
{
  class Mess{
  public:
    int FractalRank;
    int FractalNodes;
    Particle* parts_tmp;
    vector < vector <int> > Slices;
    ptrdiff_t start_x;
    ptrdiff_t length_x;
    ptrdiff_t total_memory;
    fftw_plan plan_rc;
    fftw_plan plan_cr;
    fftw_plan plan_cri;
    double* potR;
    fftw_complex* potC; 
    vector <double> green;
    MPI::Intracomm FractalWorld;
    MPI::Request send_requests[1000];
    MPI::Request recv_requests[1000];
    int send_counter;
    Mess()
    {
      cout << " made a mess " << endl;
      MPIStartup();
      FractalRank=what_is_my_rank(); 
      FractalNodes=how_many_nodes(); 
    }
    ~Mess()
    {
      FFTWFinal();
      MPIFinal();
    }
    void MPIStartup()
    {
    }
    void MPIFinal()
    {
    }
    int what_is_my_rank()
    {
    }
    int how_many_nodes()
    {
    }
    void FFTWStartup(const int& length_1,const bool& periodic)
    {
    }
    void FFTWFinal()
    {
    }
    int fftw_where(const int& i,const int& j,const int& k,const int& la,const int& lb)
    {
      return k+(j+i*la)*lb;
    }
    void calc_fftw_Slices()
    {
    }
    void How_Many_Particles_To_Send(int* countsI_out,int* countsI_in)
    {
    }
    void Send_Particles_Somewhere(int* countsI_out,int* countsI_in,
				  vector < vector <double> >& dataR_out,vector < vector <int> >& dataI_out,
				  int& how_many,double* dataR_in,int* dataI_in)
    {
    }
    void zeroR()
    {
    }
    void send_data(vector <double>& data,const int& Rank,const int& total,const bool& first,const bool& last)
    {
    }
    void recv_data(vector < vector <double> >& data,vector <int>& fromRanks,vector <int>& buffsize)
    {
    }
  };
}
#endif
