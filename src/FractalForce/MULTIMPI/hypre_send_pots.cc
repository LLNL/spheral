#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_send_pots(Fractal_Memory& mem,vector <Point*>& hypre_points,vector <double>& potH)
  {
    ofstream& FH=mem.p_file->FileHypre;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <int> counts_out(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    vector <int> counts_in(FractalNodes);
    vector <int> dataI_in;
    vector <double> dataR_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=0;
    int doubles=1;
    int FR=0;
    int off1=mem.ij_offsets[FR+1];
    const int numberB0=mem.ij_offsetsB[FractalRank];
    int numberB=numberB0;
    FH << " starting looking " << numberB << " " << mem.ij_countsB[FractalRank] << " " << potH.size() << endl;
    for(int ni=0;ni<mem.ij_countsB[FractalRank];ni++)
      {
	while(numberB >= off1)
	  {
	    FH << " looking " << FR << " " << ni << " " << numberB << " " << off1 << endl;
	    FR++;
	    assert(FR < FractalNodes);
	    off1=mem.ij_offsets[FR+1];
	  }
	//	FH << " doing " << FR << " " << ni << " " << numberB << " " << off1 << endl;
	counts_out[FR]++;
	//	dataI_out[FR].push_back(numberB-numberB0);
	dataR_out[FR].push_back(potH[ni]);
	numberB++;
      }
    //    mem.p_mess->Full_Stop();
    //    FH << " made it here in sendpotsa " << endl;
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    //    FH << " made it here in sendpotsb " << endl;
    mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    //    FH << " made it here in sendpotsc " << endl;
    assert(how_manyR == mem.ij_counts[FractalRank]);
    dataI_out.clear();
    dataR_out.clear();
    FH << " data sizeA " << how_manyR << " " << dataI_in.size() << " " << dataR_in.size() << endl;
    for(int ni=0;ni<how_manyR;ni++)
      {
	//	FH << " data sizeB " << ni << " " << dataI_in[ni] << endl;
	if(hypre_points[ni])
	  hypre_points[ni]->set_potential_point(dataR_in[ni]);
	//	hypre_points[dataI_in[ni]]->set_potential_point(dataR_in[ni]);
      }
  }
}
