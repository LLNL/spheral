#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main()
{
  MPI::Init();
  MPI::Intracomm CWorld=MPI::COMM_WORLD;
  using namespace FractalSpace;
  int procs_CWorld=CWorld.Get_size();
  int rank_CWorld=CWorld.Get_rank();
  bool master_CWorld=procs_CWorld==0;
  bool client_CWorld=!master_CWorld;
  //
  int particles;
  if(master_CWorld)
    {
      cout << " particles " << endl;
      cin >> particles;
    }



  MPI::Finalize();
  return 1;
}
namespace FractalSpace
{
  void boxx()
  {
    int splits_max=10;
    for(int split=0;split < splits_max;split++)
      {
	int dir=split % 3;
	int N=(split+1)/3;
	N=Misc::pow(2,N);
	double aN=1.0/(double)N;
	
      }
  }
}

namespace FractalSpace
{
  void counts_in_box()
  {
    for (int part=0;part < number_parts;part++)
      {

      }
  }
}
