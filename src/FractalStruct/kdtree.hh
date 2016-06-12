#ifndef _KdTree_Defined_
#define _KdTree_Defined_
namespace FractalSpace
{
  struct KdTreeNode
  {
    bool full;
    bool empty;
    int dir;
    vector <int> box;
    vector <KdTreeNode*> kids;
    list <Point*> ppoints;
    KdTreeNode();
    ~KdTreeNode();
  };
  class KdTree{
  private:
    int RANK;
    bool RANKY;
    int nnodes;
    int fullnodes;
    int spacing;
    int VOLMIN;
    double FILLFACTOR;
    KdTreeNode* rnode;
    void LoadKdTree(int corner,KdTreeNode* pnode);
    void FillBox(KdTreeNode* pnode);
    void DestroyKdTree(KdTreeNode* pnode);
    void CollectBoxes(vector <vector<int> >& SBoxes,KdTreeNode* pnode);
    void CollectPoints(vector< vector<Point*> >& SPoints,KdTreeNode* pnode);
    void CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints ,KdTreeNode* pnode);
    void Consolidate(KdTreeNode* pnode);
    void DisplayTree(KdTreeNode* pnode,int& TOT,int& NB);
    void Traverse(KdTreeNode* pnode);
  public:
    KdTree();
    ~KdTree();
    void LoadKdTree(vector <int>& BOX,vector <Point*>& pPOINTS,int spacing,int VOLMIN,double FILLFACTOR);
    void DestroyKdTree();
    void CollectBoxes(vector < vector<int> >& SBoxes);
    void CollectPoints(vector < vector<Point*> >& SPoints);
    void CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
    void Consolidate();
    void DisplayTree(int& TOT,int& NB);
    void Traverse();
  };
}
#endif
