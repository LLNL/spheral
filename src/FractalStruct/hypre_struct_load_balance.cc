#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool hypre_struct_load_balance(Fractal_Memory& mem,
				 vector<vector<int>>& SBoxes,
				 vector<vector<Point*>>& SPoints,
				 vector<int>& HRout)
  {
    const int maxLOADp=mem.hypre_max_node_load;
    const int maxLOADb=200;
    const int atLEASTp=10000;
    const int extra=500;
    const double spread=0.1;
    const int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-mem.level);
    ofstream& FHT=mem.p_file->DUMPS;
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=mem.p_mess->HypreRank;
    int HypreNodes=mem.p_mess->HypreNodes;
    const int maxtries=max(1.0,log((double)(HypreNodes)+0.1)/log(1.6));
    HRout.assign(SBoxes.size(),HypreRank);
    if(!mem.hypre_load_balance)
      return false;
    int count_points=0;
    for(auto SP : SPoints)
      count_points+=SP.size();
    vector<int>how_many;
    how_many.push_back(SPoints.size());
    how_many.push_back(count_points);
    vector<int>points_on_nodes(HypreNodes*2);
    mem.p_mess->my_AllgatherI(mem.p_mess->HypreWorld,how_many,points_on_nodes,2);
    vector<int>Boxes(HypreNodes);
    vector<int>Points(HypreNodes);
    int ni=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	Boxes[HR]=points_on_nodes[ni++];
	Points[HR]=points_on_nodes[ni++];
      }
    clean_vector(points_on_nodes);
    int total_boxes=accumulate(Boxes.begin(),Boxes.end(),0);
    int total_points=accumulate(Points.begin(),Points.end(),0);
    int average_boxes=total_boxes/HypreNodes;
    int average_points=total_points/HypreNodes;
    int mP=*max_element(Points.begin(),Points.end());
    int mB=*max_element(Boxes.begin(),Boxes.end());
    bool smallP=mP <= max(maxLOADp,atLEASTp);
    bool smallB=mB <= maxLOADb;
    FHT << " LOADY A " << mem.level << " " << mem.steps << " " << FractalRank << " " << HypreRank << " " << total_boxes << " " << total_points << " " << average_boxes << " " << average_points;
    FHT << " " << mP << " " << mB << " " << maxLOADp << " " << atLEASTp << " " << maxLOADb << " " << smallP << " " << smallB << "\n";
    if(smallP && smallB)
      return false;
    if(!smallP)
      {
	bool naughty_boy=false;
	const int SBtotal=SBoxes.size();
	for(int BOX=0;BOX<SBtotal;BOX++)
	  {
	    int over=SPoints[BOX].size()-(average_points*9)/10;;
	    if(over <= 0)
	      continue;
	    int dx=(SBoxes[BOX][1]-SBoxes[BOX][0])/spacing+1;
	    int dy=(SBoxes[BOX][3]-SBoxes[BOX][2])/spacing+1;
	    int layers=(SBoxes[BOX][5]-SBoxes[BOX][4])/spacing+1;
	    int onelayer=dx*dy;
	    int Ineed=(over-1)/onelayer+1;
	    int BZstart=layers-Ineed;
	    vector <int>bigbox=SBoxes[BOX];
	    for(int BZ=BZstart;BZ<layers;BZ++)
	      {
		SBoxes.resize(SBoxes.size()+1);
		SBoxes.back()=bigbox;
		SBoxes.back()[4]+=BZ*spacing;
		SBoxes.back()[5]=SBoxes.back()[4];
		SPoints.resize(SPoints.size()+1);
		for(int sp=BZ*onelayer;sp<(BZ+1)*onelayer;sp++)
		  SPoints.back().push_back(SPoints[BOX][sp]);
	      }
	    SPoints[BOX].resize(BZstart*onelayer);
	    SBoxes[BOX][5]=SBoxes[BOX][4]+(BZstart-1)*spacing; ///// BAAD
	    if(BZstart == 0)
	      naughty_boy=true;
	  }
	if(naughty_boy)
	  {
	    cerr << " NAUGHTY BOY A " << FractalRank << " " << mem.steps << endl;
	    int BB=0;
	    bool found_naughty_boy=false;
	    for(int B=0;B<SBoxes.size();B++)
	      {
		bool zero_box=SBoxes[B][5] < SBoxes[B][4];
		if(!zero_box && found_naughty_boy)
		  {
		    SBoxes[BB]=SBoxes[B];
		    SPoints[BB]=SPoints[B];
		  }
		if(zero_box)
		  {
		    found_naughty_boy=true;
		    cerr << " NAUGHTY BOY B " << FractalRank << " " << mem.steps << endl;
		  }
		else
		  BB++;
	      }
	    cerr << " NAUGHTY BOY C " << FractalRank << " " << mem.steps <<
	      " " << SBoxes.size() << " " << BB << endl;
	    SBoxes.resize(BB);
	    SPoints.resize(BB);
	  }
	HRout.assign(SBoxes.size(),HypreRank);
	multimap<int,deque<int>>NodesA;
	for(int HR=0;HR<HypreNodes;HR++)
	  {
	    deque<int>VHR{{HR}};
	    NodesA.insert(pair<int,deque<int>>(Points[HR],VHR));
	  }
	if((--NodesA.end())->first <= maxLOADp)
	  return false;
	int enough_spam=false;
	int tries=0;
	while((--NodesA.end())->first > maxLOADp && tries < maxtries && !enough_spam)
	  {
	    multimap<int,deque<int>>NodesB;
	    auto pstart=NodesA.begin();
	    auto pend=--NodesA.end();
	    int dd=distance(pstart,pend);
	    while(dd > 0)
	      {
		auto Va=pstart->second;
		auto Vb=pend->second;
		int aver=(Va.size()*pstart->first+Vb.size()*pend->first)/
		  (Va.size()+Vb.size());	  
		for(auto V : Vb)
		  Va.push_back(V);
		NodesB.insert(make_pair(aver,Va));
		pstart++;
		pend--;
		dd-=2;
	      }
	    if(dd == 0)
	      NodesB.insert(*pstart);
	    NodesA=NodesB;
	    tries++;
	    enough_spam=(double)(NodesB.begin()->first-(--NodesB.end())->first)/
	      (double)(NodesB.begin()->first+(--NodesB.end())->first) < spread;
	  }
	bool go_home=true;
	deque<int>mynodes;
	int aver=0;
	for(auto No : NodesA)
	  {
	    mynodes=No.second;
	    if(find(mynodes.begin(),mynodes.end(),HypreRank) != mynodes.end())
	      {
		for(auto HR : mynodes)
		  {
		    aver+=Points[HR];
		    go_home=go_home && Points[HR] <= maxLOADp;
		  }
		break;
	      }
	  }
	aver/=mynodes.size();
	if(go_home)
	  return true;
	NodesA.clear();
	int averup=(aver*102)/100;
	int averdown=(aver*98)/100;
	vector<int>can(HypreNodes,0);
	vector<int>need(HypreNodes,0);
	sort(mynodes.begin(),mynodes.end());
	for(auto HR : mynodes)
	  {
	    can[HR]=max(Points[HR]-averup,0);
	    need[HR]=max(averdown-Points[HR],0);
	  }
	vector<int>sendto(HypreNodes,0);
	for(int Hsend : mynodes)
	  {
	    if(Hsend > HypreRank)
	      break;
	    if(can[Hsend] <= 0)
	      continue;
	    for (int Hrec : mynodes)
	      {
		if(need[Hrec] <= 0 || Hrec == Hsend)
		  continue;
		int send=min(can[Hsend],need[Hrec]);
		can[Hsend]-=send;
		need[Hrec]-=send;
		if(Hsend == HypreRank)
		  sendto[Hrec]=send;
	      }
	  }
	for(auto pN=mynodes.begin();pN!=mynodes.end();pN++)
	  if(*pN == HypreRank)
	    {
	      mynodes.erase(pN);
	      break;
	    }
	sendto[HypreRank]=-1;
	int number=0;
	multimap <vector<Point*>,int,vector_comp_down>MSP;
	for(auto SP : SPoints)
	  MSP.insert(make_pair(SP,number++));
	bool fewer=false;
	for(auto M : MSP)
	  {
	    int has=M.first.size();
	    for(auto Hrec : mynodes)
	      {
		if(sendto[Hrec] > 0)
		  {
		    if(has <= sendto[Hrec]+extra)
		      {
			HRout[M.second]=Hrec;
			sendto[Hrec]-=has;
			break;
		      }
		  }
		else
		  fewer=true;
	      }
	    if(!fewer)
	      continue;
	    auto pN=mynodes.begin();
	    while(pN != mynodes.end())
	      {
		if(sendto[*pN] <= 0)
		  pN=mynodes.erase(pN);
		else
		  pN++;
	      }
	  }
      }
    if(!smallB)
      move_small_boxes(mem,Boxes,SPoints,HRout);
    return true;
  }
}
