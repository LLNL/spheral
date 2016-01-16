#ifndef _Octree_Defined_
#define _Octree_Defined_
namespace FractalSpace
{
  struct OcTreeNode
  {
    bool full;
    int* box;
    OcTreeNode** kids;
    list <Point*> ppoints;
    OcTreeNode();
    ~OcTreeNode();
  };
  class OcTree{
  private:
    int VolumeFactor;
    OcTreeNode* rnode;
    void LoadOcTree(int corner,OcTreeNode* pnode);
    void DestroyOcTree(OcTreeNode* pnode);
    void CollectBoxes(vector <vector<int> >& SSBoxes,OcTreeNode* pnode);
    void CollectPoints(vector< vector<Point*> >& SSPoints,OcTreeNode* pnode);
    void CollectBoxesPoints(vector < vector<int> >& SSBoxes,vector < vector<Point*> >& SSPoints ,OcTreeNode* pnode);
    void Traverse(OcTreeNode* pnode);
  public:
    OcTree();
    ~OcTree();
    void LoadOcTree(int* BOX,vector <Point*>& pPOINTS,int spacing);
    void DestroyOcTree();
    void CollectBoxes(vector < vector<int> >& SSBoxes);
    void CollectPoints(vector < vector<Point*> >& SSPoints);
    void CollectBoxesPoints(vector < vector<int> >& SSBoxes,vector < vector<Point*> >& SSPoints);
    void Traverse();
  };
}
#endif
