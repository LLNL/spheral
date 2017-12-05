#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void node_groups_struct(Fractal_Memory& mem,vector <int>& counts)
  {
    static int COUNTER=0;
    ofstream& FHT=mem.p_file->DUMPS;
    int FractalNodes=mem.p_mess->FractalNodes;
    int FractalRank=mem.p_mess->FractalRank;
    mem.p_mess->Full_Stop_Do_Not_Argue();
    bool allfull=*min_element(counts.begin(),counts.end()) > 0;
    FHT << " FULL " << mem.level << " " << mem.steps << " " << allfull << "\n";
    vector <int> head_number;
    if(allfull)
      head_number.resize(FractalNodes,0);
    else
      {
	vector <int> counts_in(FractalNodes);
	vector <int> counts_out(FractalNodes,0);
	vector <int> dataI_in;
	vector <double> dataR_in;
	vector < vector <int> > dataI_out(FractalNodes);
	vector < vector <double> > dataR_out(FractalNodes);
	bool tryit=mem.p_mess->IAmAHypreNode;
	for(int FR : mem.Touchy)
	  if(tryit && FR > FractalRank)
	    dataI_in.push_back(FR);
	int ss=dataI_in.size();
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    if(tryit)
	    // if(tryit && counts[FR] > 0)
	      {
		dataI_out[FR]=dataI_in;
		counts_out[FR]=ss;
	      }
	  }
	clean_vector(dataI_in);
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=1;
	int doubles=0;
	mem.p_mess->Send_Data_Some_How(8,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	clean_vector(counts_out);
	dataI_out.clear();
	dataR_out.clear();
	vector <int> list_pair_1;
	vector <int> list_pair_2;
	int counterI=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    head_number.push_back(FR);
	    for(int c=0;c<counts_in[FR];c++)
	      {
		list_pair_1.push_back(FR);
		list_pair_2.push_back(dataI_in[counterI]);
		counterI++;
	      }
	  }
	clean_vector(dataI_in);
	clean_vector(dataR_in);
	for(int i=0; i < counterI; ++i)
	  {
	    int j=list_pair_1[i];
	    while(head_number[j] != j)
	      j=head_number[j];
	    int k=list_pair_2[i];
	    while(head_number[k] != k)
	      k=head_number[k];
	    if(j != k) 
	      head_number[j]=k;
	  }
	clean_vector(list_pair_1);
	clean_vector(list_pair_2);
	for (int j=0;j < FractalNodes; ++j)
	  {
	    while(head_number[j] != head_number[head_number[j]])
	      head_number[j]=head_number[head_number[j]];
	  }
      }
    clean_vector(mem.p_mess->Hranks);
    int myhead=head_number[FractalRank];
    mem.p_mess->IHranks.resize(FractalNodes);
    for(int FR=0;FR<FractalNodes;FR++)
      if(head_number[FR] == myhead)
	{
	  mem.p_mess->Hranks.push_back(FR);
	  mem.p_mess->IHranks[FR]=mem.p_mess->Hranks.size()-1;
	}
    if(!mem.p_mess->IAmAHypreNode)
      {
	clean_vector(mem.p_mess->Hranks);
	mem.p_mess->IHranks.assign(FractalNodes,-1);
      }
    mem.p_mess->HypreNodes=mem.p_mess->Hranks.size();
    mem.p_mess->node_lists.clear();
    clean_vector(mem.p_mess->freenodes);
    if(mem.hypre_load_balance)
      {
	{
	  map<int,int>headmap;
	  for(int FR=0;FR<FractalNodes;FR++)
	    if(counts[FR] > 0 && head_number[FR] == FR)
	      headmap.insert(make_pair(FR,headmap.size()));
	  mem.p_mess->node_lists.resize(headmap.size());
	  for(int FR=0;FR<FractalNodes;FR++)
	    {
	      if(counts[FR] == 0)
		mem.p_mess->freenodes.push_back(FR);
	      else
		{
		  auto it=headmap.find(head_number[FR]);
		  assert(it != headmap.end());
		  mem.p_mess->node_lists[it->second].push_back(FR);
		}
	    }
	}
	clean_vector(head_number);
	int nla=0;
	int nlb=0;
	for(auto &nl : mem.p_mess->node_lists)
	  {
	    if(nl.size() > 1)
	      {
		if(nla < nlb)
		  mem.p_mess->node_lists[nla]=mem.p_mess->node_lists[nlb];
		nla++;
	      }
	    nlb++;
	  }
	mem.p_mess->node_lists.resize(nla);
	int nln=0;
	for(auto nl : mem.p_mess->node_lists)
	  {
	    int nlr=0;
	    for(auto FR : nl)
	      {
		FHT << " NODELIST A " << mem.steps << " " << mem.level << " " << FractalRank << " " << nln;
		FHT << " " << nlr++ << " " << FR << "\n";
	      }
	    nln++;
	  }
	int fnl=0;
	for(auto FN : mem.p_mess->freenodes)
	  FHT << " FREE A " << mem.steps << " " << mem.level << " " << FractalRank << " " << fnl++ << " " << FN << "\n";
      }
    COUNTER++;
  }
}
