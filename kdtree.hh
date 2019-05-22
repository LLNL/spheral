#ifndef _KdTree_Defined_
#define _KdTree_Defined_
namespace FractalSpace
{
  struct KdTreeNode
  {
    bool full;
    bool empty;
    int dir;
    std::vector <int> box;
    std::vector <KdTreeNode*> kids;
    std::list <Point*> ppoints;
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
    void CollectBoxes(std::vector <std::vector<int> >& SBoxes,KdTreeNode* pnode);
    void CollectPoints(std::vector< std::vector<Point*> >& SPoints,KdTreeNode* pnode);
    void CollectBoxesPoints(std::vector < std::vector<int> >& SBoxes,std::vector < std::vector<Point*> >& SPoints ,KdTreeNode* pnode);
    void Consolidate(KdTreeNode* pnode);
    void DisplayTree(KdTreeNode* pnode,int& TOT,int& NB);
    void Traverse(KdTreeNode* pnode);
  public:
    KdTree();
    ~KdTree();
    void LoadKdTree(std::vector <int>& BOX,std::vector <Point*>& pPOINTS,int spacing,int VOLMIN,double FILLFACTOR);
    void DestroyKdTree();
    void CollectBoxes(std::vector < std::vector<int> >& SBoxes);
    void CollectPoints(std::vector < std::vector<Point*> >& SPoints);
    void CollectBoxesPoints(std::vector < std::vector<int> >& SBoxes,std::vector < std::vector<Point*> >& SPoints);
    void Consolidate();
    void DisplayTree(int& TOT,int& NB);
    void Traverse();
  };
}
#endif
