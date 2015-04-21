#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void node_groups(Fractal_Memory& mem)
  {
    ofstream& FF=mem.p_mess->p_file->DUMPS;
    int FractalNodes=mem.p_mess->FractalNodes;
    int FractalRank=mem.p_mess->FractalRank;
    //    cerr << " NODENODE " << FractalRank << " " << mem.Touchy.size() << "\n";
    vector <int> counts_in(FractalNodes);
    vector <int> counts_out(FractalNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);
    bool tryit=mem.ij_counts[FractalRank] > 0;
    for(int T=0;T<mem.Touchy.size();T++)
      {
	//	FF << " TOUCHY0 " << mem.Touchy[T] << "\n";
	if(tryit && mem.Touchy[T] > FractalRank)
	  {
	    dataI_in.push_back(mem.Touchy[T]);
	    //	    FF << " TOUCHYA " << dataI_in.back() << "\n";
	  }
      }
    int ss=dataI_in.size();
    for(int FR=0;FR<FractalNodes;FR++)
      {
	if(tryit && mem.ij_counts[FR] > 0)
	  {
	    dataI_out[FR]=dataI_in;
	    counts_out[FR]=ss;
	  }
      }
    dataI_in.clear();
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=0;
    mem.p_mess->Send_Data_Some_How(8,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    vector <int> list_pair_1;
    vector <int> list_pair_2;
    vector <int> head_number;
    int counterI=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	head_number.push_back(FR);
	for(int c=0;c<counts_in[FR];c++)
	  {
	    list_pair_1.push_back(FR);
	    list_pair_2.push_back(dataI_in[counterI]);
	    //	    FF << " TOUCHYB " << list_pair_1.back() << " " << list_pair_2.back() << "\n";
	    counterI++;
	  }
      }
    dataI_in.clear();
    dataR_in.clear();
    for(int i=0; i < counterI; ++i)
      {
	int j=list_pair_1[i];
	while(head_number[j]  != j)
	  j=head_number[j];
	int k=list_pair_2[i];
	while(head_number[k] != k)
	  k=head_number[k];
	if(j != k) 
	  head_number[j]=k;
      }
    list_pair_1.clear();
    list_pair_2.clear();
    for (int j=0;j < FractalNodes; ++j)
      {
	while(head_number[j] != head_number[head_number[j]])
	  head_number[j]=head_number[head_number[j]];
      }
    mem.Touchy.clear();
    int myhead=head_number[FractalRank];
    for(int FR=0;FR<FractalNodes;FR++)
      {
	//	FF << " NODESA " << FR << " " << head_number[FR] << " " << myhead << "\n";
	if(head_number[FR] != myhead)
	  mem.ij_counts[FR]=0;
	else
	  {
	    mem.Touchy.push_back(FR);
	    //	    FF << " NODESB " << FR << " " << mem.ij_counts[FR] << "\n";
	  }
	//	FF << " NODESC " << FR << " " << mem.ij_counts[FR] << "\n";
      }
  }
}
