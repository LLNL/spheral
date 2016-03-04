#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_boxes(vector <vector <Point*> >hypre_points,int spacing,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
//     int RANK=-1;
//     MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
//     bool Ranky= RANK == 21;
    cerr << " Enter BOXES A " << RANK << " " << hypre_points.size() << endl;
    int MAXY=Misc::pow(2,29);
    int MINY=-Misc::pow(2,29);
    int ni=0;
    MPI_Barrier(MPI_COMM_WORLD);
    vector <int>BOX(6);
    for(vector <Point*>& hp : hypre_points)
      {
	vector <int> BOX(6);
	vector <int> pos(3);
	for(int B : {0,1,2,3,4,5})
	  BOX[B]=B%2 == 0 ? MAXY:MINY;
	for(Point* &p : hp)
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
	OcTree* pHypTree=new OcTree();
	pHypTree->LoadOcTree(BOX,hypre_points[ni],spacing);
	int TotalPoints=0;
	int TotalBoxes=0;
	pHypTree->DisplayTree(TotalPoints,TotalBoxes);
// 	cerr << " BOXTotal " << RANK << " " << ni << " " << hypre_points[ni].size() << " " << TotalPoints << " " << TotalBoxes << endl;
	pHypTree->CollectBoxesPoints(SBoxes,SPoints);
	delete pHypTree;
	ni++;
      }
    for(vector<int>& SB : SBoxes)
      {
	for(int B : {1,3,5})
	  SB[B]--;
	Misc::times(SB,spacing);
      }
  }
}
