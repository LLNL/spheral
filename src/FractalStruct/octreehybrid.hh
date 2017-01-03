#ifndef _Octree_Defined_
#define _Octree_Defined_
namespace FractalSpace
{
  struct OcTreeNode
  {
    bool full;
    bool empty;
    int dir;
    vector <int> box;
    vector <OcTreeNode*> kids;
    list <Point*> ppoints;
    OcTreeNode();
    ~OcTreeNode();
  };
  class OcTree{
  private:
    int RANK;
    int nnodes;
    int fullnodes;
    int spacing;
    int VOLMIN;
    double FILLFACTOR;
    OcTreeNode* rnode;
    void LoadOcTree(int corner,OcTreeNode* pnode);
    void DestroyOcTree(OcTreeNode* pnode);
    void CollectBoxes(vector <vector<int> >& SBoxes,OcTreeNode* pnode);
    void CollectPoints(vector< vector<Point*> >& SPoints,OcTreeNode* pnode);
    void CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints ,OcTreeNode* pnode);
    void Consolidate(OcTreeNode* pnode);
    void DisplayTree(OcTreeNode* pnode,int& TOT,int& NB);
    void Traverse(OcTreeNode* pnode);
  public:
    OcTree();
    ~OcTree();
    void LoadOcTree(vector <int>& BOX,vector <Point*>& pPOINTS,int spacing,int VOLMIN,double FILLFACTOR);
    void DestroyOcTree();
    void CollectBoxes(vector < vector<int> >& SBoxes);
    void CollectPoints(vector < vector<Point*> >& SPoints);
    void CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
    void Consolidate();
    void DisplayTree(int& TOT,int& NB);
    void Traverse();
  };
}
#endif
