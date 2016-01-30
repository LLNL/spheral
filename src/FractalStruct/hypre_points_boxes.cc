#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_boxes(vector <vector <Point*> >hypre_points,int spacing,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    int BOX[6];
    vector <int>pos(3);
    for(int ni=0;ni<hypre_points.size();ni++)
      {
	for(int B=0;B<6;B++)
	  BOX[B]=B%2 == 0 ? INT_MAX:INT_MIN;
	for(vector<Point*>::const_iterator point_itr=hypre_points[ni].begin();
	    point_itr !=hypre_points[ni].end();++point_itr)
	  {
	    (*point_itr)->get_pos_point(pos);
	    for(int B=0;B<6;B+=2)
	      {
		BOX[B]=min(BOX[B],pos[B/2]);
		BOX[B+1]=max(BOX[B+1],pos[B/2]);
	      }
	  }
	for(int B=1;B<6;B+=2)
	  BOX[B]+=spacing;
	OcTree* pHypTree=new OcTree();
	pHypTree->LoadOcTree(BOX,hypre_points[ni],spacing);
	pHypTree->CollectBoxesPoints(SBoxes,SPoints);
	pHypTree->DestroyOcTree();
	delete pHypTree;
      }
  }
}
