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
      cout << " cleaned up a mess " << endl;
      FFTWFinal();
      MPIFinal();
    }
    void MPIStartup()
    {
      if(!MPI::Is_initialized())
	MPI::Init();
      FractalWorld=MPI::COMM_WORLD;
    }
    void MPIFinal()
    {
      //      if(!MPI::Is_finalized())
      MPI::Finalize();
    }
    int what_is_my_rank()
    {
      return FractalWorld.Get_rank();
    }
    int how_many_nodes()
    {
      return FractalWorld.Get_size();
    }
    void FFTWStartup(const int& length_1,const bool& periodic)
    {
      fftw_mpi_init();
      if(periodic)
	{
	  const int length_c=(length_1+2)/2;
	  total_memory=fftw_mpi_local_size_3d(length_1,length_1,length_c,FractalWorld,&length_x,&start_x);
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(length_1,length_1,length_1,potR,potC,FractalWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(length_1,length_1,length_1,potC,potR,FractalWorld,FFTW_MEASURE);
	  plan_cri=fftw_mpi_plan_dft_c2r_3d(length_1,length_1,length_1,potC,potR,FractalWorld,FFTW_ESTIMATE);
	}
      else
	{
	  const int length_11=length_1+1;
	  const int length_2=2*length_1;
	  const double g_c=pow(static_cast<double>(length_1),-5)/8.0;
	  total_memory=fftw_mpi_local_size_3d(length_2,length_2,length_11,FractalWorld,&length_x,&start_x);
	  green.resize(static_cast<int>(length_x)*length_11*length_11);
	  size_t sizeR=sizeof(double);
	  size_t sizeC=sizeof(fftw_complex);
	  double* greenR;
	  fftw_complex* greenC;
	  greenR=(double*) fftw_malloc(sizeR*2*total_memory);
	  greenC=(fftw_complex*) fftw_malloc(sizeC*total_memory);
	  //	  greenR=fftw_alloc_real(2*total_memory);   ***********
	  //	  greenC=fftw_alloc_complex(total_memory);  ******************
	  fftw_plan plan_green_rc=fftw_mpi_plan_dft_r2c_3d(length_2,length_2,length_2,greenR,greenC,FractalWorld,FFTW_ESTIMATE);
	  for(int n_x=0;n_x < length_x;++n_x)
	    {
	      int n_xa=n_x+start_x;
	      int n_x_inv=length_2-n_xa;
	      n_xa=min(n_xa,n_x_inv);
	      double x_2=static_cast<double>(n_xa*n_xa);
	      for(int n_y=0;n_y < length_11;++n_y)
		{
		  int n_y_inv=n_y;
		  if(n_y > 0) n_y_inv=length_2-n_y;
		  double y_2=static_cast<double>(n_y*n_y);
		  for( int n_z=0;n_z<length_11;++n_z)
		    {
		      int n_z_inv=n_z;
		      if(n_z > 0) n_z_inv=length_2-n_z;       
		      double z_2=static_cast<double>(n_z*n_z)+0.25;
		      double r_2=z_2+y_2+x_2;
		      double dead_parrot=-g_c/sqrt(r_2);
		      greenR[fftw_where(n_x,n_y,n_z,length_2,length_2)]=dead_parrot;
		      greenR[fftw_where(n_x,n_y_inv,n_z,length_2,length_2)]=dead_parrot;
		      greenR[fftw_where(n_x,n_y,n_z_inv,length_2,length_2)]=dead_parrot;
		      greenR[fftw_where(n_x,n_y_inv,n_z_inv,length_2,length_2)]=dead_parrot;
		    }	      
		}
	    }
	  fftw_execute(plan_green_rc);
	  for( int p_x=0;p_x<length_x;++p_x)
	    {
	      for( int p_y=0;p_y<length_11;++p_y)
		{
		  for( int p_z=0;p_z<length_11;++p_z)
		    {
		      green[fftw_where(p_x,p_y,p_z,length_11,length_11)]=greenC[fftw_where(p_x,p_y,p_z,length_2,length_11)][0];
		    }
		}
	    }
	  fftw_free(greenR);
	  fftw_free(greenC);
	  greenR=0;
	  greenC=0;
	  fftw_destroy_plan(plan_green_rc);
	  //
	  plan_rc=fftw_mpi_plan_dft_r2c_3d(length_2,length_2,length_2,potR,potC,FractalWorld,FFTW_MEASURE);
	  plan_cr=fftw_mpi_plan_dft_c2r_3d(length_2,length_2,length_2,potC,potR,FractalWorld,FFTW_MEASURE);
	}
    }
    void FFTWFinal()
    {
      fftw_free(potR);
      fftw_free(potC);
      fftw_destroy_plan(plan_rc);
      fftw_destroy_plan(plan_cr);
      fftw_mpi_cleanup();
    }
    int fftw_where(const int& i,const int& j,const int& k,const int& la,const int& lb)
    {
      return k+(j+i*la)*lb;
    }
    void calc_fftw_Slices()
    {
      vector <int>paramsend(2);
      vector <int>paramrecv(2*FractalNodes);
      paramsend[0]=start_x;
      paramsend[1]=start_x+length_x-1;
      FractalWorld.Allgather(&paramsend,2,MPI::INT,&paramrecv,2,MPI::INT);
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  Slices[FR].resize(2);
	  Slices[FR][0]=paramrecv[2*FR];
	  Slices[FR][1]=paramrecv[2*FR+1];
	}
    }
    void How_Many_Particles_To_Send(int* countsI_out,int* countsI_in)
    {
      for(int FR=0;FR<FractalNodes;FR++)
	FractalWorld.Gather(&countsI_out[FR],1,MPI::INT,countsI_in,1,MPI::INT,FR);
    }
    void Send_Particles_Somewhere(int* countsI_out,int* countsI_in,
				  vector < vector <double> >& dataR_out,vector < vector <int> >& dataI_out,
				  int& how_many,double* dataR_in,int* dataI_in)
    {
      int* displsI= new int[FractalNodes];
      int* displsR= new int[FractalNodes];
      int* countsR_in= new int[FractalNodes];
      int* countsR_out= new int[FractalNodes];
      displsI[0]=0;
      displsR[0]=0;
      for(int FR=1;FR<FractalNodes;FR++)
	{
	  displsI[FR]=displsI[FR-1]+countsI_in[FR-1];
	  displsR[FR]=4*displsI[FR];
	  countsR_out[FR]=countsI_out[FR]*4;
	  countsR_in[FR]=countsI_in[FR]*4;
	}
      how_many=displsI[FractalNodes-1]+countsI_in[FractalNodes-1];
      dataR_in=new double[how_many];
      dataI_in=new int[how_many];
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  FractalWorld.Gatherv(&dataR_out[FR],countsR_out[FR],MPI::DOUBLE,dataR_in,countsR_in,displsR,MPI::DOUBLE,FR);  
	  FractalWorld.Gatherv(&dataI_out[FR],countsI_out[FR],MPI::INT,dataI_in,countsI_in,displsI,MPI::INT,FR);  
	}

      delete [] displsI;
      delete [] displsR;
      delete [] countsR_in;
      delete [] countsR_out;
      displsI=0;
      displsR=0;
      countsR_in=0;
      countsR_out=0;
      dataI_out.clear();
      dataR_out.clear();
    }
    void zeroR()
    {
      for(ptrdiff_t ni=0;ni<total_memory;ni++)
	potR[ni]=0.0;
    }

    void send_data(vector <double>& data,const int& Rank,const int& total,const bool& first,const bool& last)
    {
      if(first)
	send_counter=0;
      send_requests[send_counter] = FractalWorld.Isend(&data,total,MPI::DOUBLE,Rank,0);
      send_counter++;
      if(last)
	{
	  MPI::Request black_knight;
	  black_knight.Waitall(send_counter,send_requests);
	  send_counter=0;
	}
    }
    void recv_data(vector < vector <double> >& data,vector <int>& fromBoxes,vector <int>& buffsize)
    {
      int SBNodes=fromBoxes.size();
      for(int SB=0;SB<SBNodes;SB++)
	recv_requests[SB] = FractalWorld.Irecv(&data[SB],buffsize[SB],MPI::DOUBLE,fromBoxes[SB],0);
      MPI::Request black_knight;
      black_knight.Waitall(SBNodes,recv_requests);
    }
  };
}
#endif
