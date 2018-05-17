#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_boxes(Fractal_Memory& mem,vector <vector <Point*> >& hypre_points,int spacing,
			  int VOLMIN,double FILLFACTOR,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    static int _COUNTER=0;
    // ofstream& FHT=mem.p_file->DUMPS;
    const int SBstart=SBoxes.size();
    int RANK=-1;
    MPI_Comm_rank(Fractal_Memory::FRACTAL_UNIVERSE,&RANK);
    // bool RANKY=RANK==21;
    const int MAXY=Misc::pow(2,29);
    const int MINY=-MAXY;
    for(auto hp : hypre_points)
      {
	vector <int> BOX({{MAXY,MINY,MAXY,MINY,MAXY,MINY}});
	vector <int> pos(3);
	for(auto p : hp)
	  {
	    p->get_pos_point(pos);
	    for(int B : {0,2,4})
	      {
		BOX[B]=min(BOX[B],pos[B/2]);
		BOX[B+1]=max(BOX[B+1],pos[B/2]);
	      }
	  }
	Misc::divide(BOX,spacing);
	for(int B : {1,3,5})
	  BOX[B]++;
	KdTree* pHypTree=new KdTree();
	pHypTree->LoadKdTree(BOX,hp,spacing,VOLMIN,FILLFACTOR);
	int TotalPoints=0;
	int TotalBoxes=0;
	pHypTree->DisplayTree(TotalPoints,TotalBoxes);
	pHypTree->CollectBoxesPoints(SBoxes,SPoints);
	delete pHypTree;
      }
    for(int SBcount=SBstart;SBcount<SBoxes.size();SBcount++)
      {
	for(int B : {1,3,5})
	  SBoxes[SBcount][B]--;
	Misc::times(SBoxes[SBcount],spacing);
      }
    _COUNTER++;
    if(VOLMIN > 1 || FILLFACTOR < 1.0)
      any_overlaps(mem,spacing,VOLMIN,FILLFACTOR,SBoxes,SPoints);
    auto SP=SPoints.begin();
    for(auto SB : SBoxes)
      {
	int vol=(SB[1]-SB[0])/spacing+1;
	vol*=(SB[3]-SB[2])/spacing+1;
	vol*=(SB[5]-SB[4])/spacing+1;
	int volS=(*SP).size();
	if(volS != vol || vol == 0 || volS == 0)
	  {
	    cerr << " BBBBAAADD BOXES " << RANK << " " << vol << " " << volS << endl;
	    cerr << SB[0] << " " << SB[1] << " " << SB[2] << " " << SB[3] << " " << SB[4] << " " << SB[5] << endl;
	  }
	SP++;
      }
  }
}
