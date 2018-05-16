#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  KdTreeNode::KdTreeNode():
    full(false),
    empty(false)
  {}
  KdTreeNode::~KdTreeNode()
  {}
  KdTree::KdTree():
    nnodes(0),
    fullnodes(0),
    rnode(NULL),
    spacing(1),
    VOLMIN(-1),
    FILLFACTOR(2.0)
  {
    RANK=-1;
    MPI_Comm_rank(Fractal_Memory::FRACTAL_UNIVERSE,&RANK);
    RANKY=RANK==-21;
  }
  KdTree::~KdTree()
  {
    DestroyKdTree();
  }
  void KdTree::LoadKdTree(vector<int>& BOX,vector <Point*>& pPOINTS,int space=1,int VMIN=-1,double FFACTOR=2.0)
  {
    nnodes++;
    spacing=space;
    VOLMIN=VMIN;
    FILLFACTOR=FFACTOR;
    if(rnode == NULL)
      {
	try
	  {
	    rnode=new KdTreeNode;
	  }
	catch(bad_alloc& ba)
	  {
	    cerr << " badd root node " << ba.what() << " " << RANK << endl;
	    exit(0);
	  }
      }
    rnode->dir=-1;
    rnode->kids.assign(2,NULL);
    rnode->box=BOX;
    rnode->ppoints.assign(pPOINTS.begin(),pPOINTS.end());
    rnode->empty=rnode->ppoints.empty();
    if(rnode->empty)
      {
	rnode->full=false;
	return;
      }
    int VOL=(BOX[1]-BOX[0])*(BOX[3]-BOX[2])*(BOX[5]-BOX[4]);
    rnode->full=rnode->ppoints.size() == VOL;
    const double V8=0;
    if(VOL > 1 && !rnode->full && (VOLMIN > 1 || FILLFACTOR < 1.0))
      {
	double ff=(double)(rnode->ppoints.size())/(double)(VOL);
	bool veryFULL= ff >= FILLFACTOR;
	bool smallVOL=(VOL <= VOLMIN && ff > 0.249) || VOL <= V8;
	if(smallVOL || veryFULL  )
	  FillBox(rnode);
      }
    rnode->full=rnode->ppoints.size() == VOL;
    if(rnode->full)
      {
	fullnodes++;
	return;
      }
    for(int ni : {0,1})
      if(rnode->kids[ni] == NULL)
	{
	  try
	    {
	      rnode->kids[ni]= new KdTreeNode;
	    }
	  catch(bad_alloc& ba)
	    {
	      cerr << " badd root kid " << ba.what() << " " << ni << " " << RANK << endl;
	      exit(0);
	    }
	}
    for(int corner : {0,1}) 
      {
	rnode->kids[corner]->dir=2;
	LoadKdTree(corner,rnode);
      }
    nnodes=0;
    fullnodes=0;
  }
  void KdTree::LoadKdTree(int cornera,KdTreeNode* pnode)
  {
    static const int MINX=Misc::pow(2,29);
    static const int MAXX=-MINX;
    nnodes++;
    KdTreeNode* knode=pnode->kids[cornera];
    knode->kids.assign(2,NULL);
    knode->box=pnode->box;
    int dir=knode->dir;
    int DIR2=2*dir;
    knode->box[cornera % 2 == 0 ? DIR2+1:DIR2]=(pnode->box[DIR2]+pnode->box[DIR2+1])/2;
    knode->ppoints.clear();
    vector <int>pos(3);
    vector<int>KBOXA={MINX,MINX,MINX};
    vector<int>KBOXB={MAXX,MAXX,MAXX};
    auto itp=pnode->ppoints.begin();
    auto itpe=pnode->ppoints.end();
    while(itp!=itpe)
      {
	(*itp)->get_pos_point(pos);
   	Misc::divide(pos,spacing);
	if(pos[0] < knode->box[0] || pos[1] < knode->box[2] || pos[2] < knode->box[4] ||
	   pos[0] >= knode->box[1] || pos[1] >= knode->box[3] || pos[2] >= knode->box[5])
	  itp++;
	else
	  {
	    knode->ppoints.push_back(*itp);
	    itp=pnode->ppoints.erase(itp);
	    for(int ni : {0,1,2}) 
	      {
		KBOXA[ni]=min(KBOXA[ni],pos[ni]);
		KBOXB[ni]=max(KBOXB[ni],pos[ni]);
	      }
	  }
      }
    int ni2=0;
    for(int ni : {0,1,2}) 
      {
	knode->box[ni2]=KBOXA[ni];
	knode->box[ni2+1]=KBOXB[ni]+1;
	ni2+=2;
      }
    int vol=(knode->box[1]-knode->box[0]);
    vol*=(knode->box[3]-knode->box[2]);
    vol*=(knode->box[5]-knode->box[4]);
    knode->empty=vol == 0 || knode->ppoints.empty();
    if(knode->empty)
      {
	knode->full=false;
	return;
      }
    knode->full=vol != 0 && knode->ppoints.size() == vol;
    const double V8=0;
    double ff=(double)(knode->ppoints.size())/(double)(vol);
    bool veryFULL= ff >= FILLFACTOR;
    bool smallVOL=(vol <= VOLMIN && ff > 0.249) || vol <= V8;
    if((VOLMIN > 1 || FILLFACTOR < 1.0) && !knode->full && vol > 1 && (smallVOL || veryFULL))
      FillBox(knode);
    knode->full=vol != 0 && knode->ppoints.size() == vol;
    if(knode->full)
      {
	fullnodes++;
      }
    if(knode->full)
      return;
    for(int ni : {0,1}) 
      if(knode->kids[ni] == NULL)
	{
	  try
	    {
	      knode->kids[ni]= new KdTreeNode;
	    }
	  catch(bad_alloc& ba)
	    {
	      cerr << " badd node kid " << ba.what() << " " << ni << " " << RANK << endl;
	      exit(0);
	    }
	}
    for(int cornerb : {0,1}) 
      {
	knode->kids[cornerb]->dir=(knode->dir+1)%3;
	LoadKdTree(cornerb,knode);
      }
  }
  void KdTree::FillBox(KdTreeNode* pnode)
  {
    int vol=(pnode->box[1]-pnode->box[0]);
    vol*=(pnode->box[3]-pnode->box[2]);
    vol*=(pnode->box[5]-pnode->box[4]);
    if(RANKY)
      cerr << " FILLA " << vol << " " << VOLMIN << " " << FILLFACTOR << " " << pnode->ppoints.size() << "\n";
    if(vol == pnode->ppoints.size())
      return;
    vector<int>pos(3);
    std::map<array<int,3>,Point*,point_comp2> boxP;
    for(auto pp : pnode->ppoints)
      {
	auto ret=boxP.insert(make_pair(pp->get_pos_point_a(),pp));
	assert(ret.second);
      }
    if(RANKY)
      cerr << " FILLB " << vol << " " << VOLMIN << " " << FILLFACTOR << " " << pnode->ppoints.size() << "\n";
    pnode->ppoints.clear();
    Point* pFAKE=0;
    for(int nz=pnode->box[4];nz<pnode->box[5];nz++)
      {
	int pz=nz*spacing;
	for(int ny=pnode->box[2];ny<pnode->box[3];ny++)
	  {
	    int py=ny*spacing;
	    for(int nx=pnode->box[0];nx<pnode->box[1];nx++)
	      {
		int px=nx*spacing;
		array<int,3>posa3{{px,py,pz}};
		if(boxP.find(posa3) == boxP.end())
		   boxP[posa3]=pFAKE;
	      }
	  }
      }
    for(auto mP : boxP)
      pnode->ppoints.push_back(mP.second);
    if(RANKY)
      cerr << " FILLC " << vol << " " << VOLMIN << " " << FILLFACTOR << " " << pnode->ppoints.size() << "\n";
  }
  void KdTree::DestroyKdTree()
  {
    DestroyKdTree(rnode);
  }
  void KdTree::DestroyKdTree(KdTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    DestroyKdTree(pnode->kids[0]);
    DestroyKdTree(pnode->kids[1]);
    delete pnode;
  }
  void KdTree::CollectBoxes(vector < vector<int> >& SBoxes)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=rnode->box;
	return;
      }
    CollectBoxes(SBoxes,rnode->kids[0]);
    CollectBoxes(SBoxes,rnode->kids[1]);
  }
  void KdTree::CollectBoxes(vector< vector<int> >& SBoxes,KdTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=pnode->box;
	return;
      }
    CollectBoxes(SBoxes,pnode->kids[0]);
    CollectBoxes(SBoxes,pnode->kids[1]);
  }
  void KdTree::CollectPoints(vector < vector<Point*> >& SPoints)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    CollectPoints(SPoints,rnode->kids[0]);
    CollectPoints(SPoints,rnode->kids[1]);
  }
  void KdTree::CollectPoints(vector< vector<Point*> >& SPoints,KdTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    CollectPoints(SPoints,pnode->kids[0]);
    CollectPoints(SPoints,pnode->kids[1]);
  }
  void KdTree::CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=rnode->box;
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    CollectBoxesPoints(SBoxes,SPoints,rnode->kids[0]);
    CollectBoxesPoints(SBoxes,SPoints,rnode->kids[1]);
  }
  void KdTree::CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints,KdTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=pnode->box;
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    CollectBoxesPoints(SBoxes,SPoints,pnode->kids[0]);
    CollectBoxesPoints(SBoxes,SPoints,pnode->kids[1]);
  }
  void KdTree::DisplayTree(int& TOT,int& NB)
  {
    if(rnode != NULL)
      {
	int vol=(rnode->box[1]-rnode->box[0]);
	vol*=(rnode->box[3]-rnode->box[2]);
	vol*=(rnode->box[5]-rnode->box[4]);
	if(rnode->full)
	  {
	    assert(rnode->ppoints.size() == vol);
	    assert(vol != 0);
	    TOT+=rnode->ppoints.size();
	    NB++;
	    return;
	  }
	else
	  assert(rnode->ppoints.size() < vol);
	DisplayTree(rnode,TOT,NB);
      }
  }
  void KdTree::DisplayTree(KdTreeNode* pnode,int& TOT,int& NB)
  {
    if(pnode != NULL)
      {
	int vol=(pnode->box[1]-pnode->box[0]);
	vol*=(pnode->box[3]-pnode->box[2]);
	vol*=(pnode->box[5]-pnode->box[4]);
	if(pnode->full)
	  {
	    assert(pnode->ppoints.size() == vol);
	    assert(vol != 0);
	    TOT+=pnode->ppoints.size();
	    NB++;
	    return;
	  }
	else
	  if(vol > 0)
	    assert(pnode->ppoints.size() < vol);
	DisplayTree(pnode->kids[0],TOT,NB);
	DisplayTree(pnode->kids[1],TOT,NB);
      }
  }
  void KdTree::Traverse()
  {
    if(rnode != NULL)
      Traverse(rnode);
  }
  void KdTree::Traverse(KdTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    Traverse(pnode->kids[0]);
    Traverse(pnode->kids[1]);
  }
}
