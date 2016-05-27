#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  OcTreeNode::OcTreeNode():
    full(false),
    empty(false)
  {}
  OcTreeNode::~OcTreeNode()
  {}
  OcTree::OcTree():
    nnodes(0),
    fullnodes(0),
    rnode(NULL),
    spacing(1),
    VOLMIN(-1),
    FILLFACTOR(2.0)
  {
    RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
  }
  OcTree::~OcTree()
  {
    DestroyOcTree();
  }
  void OcTree::LoadOcTree(vector<int>& BOX,vector <Point*>& pPOINTS,int space=1,int VMIN=-1,double FFACTOR=2.0)
  {
    nnodes++;
    spacing=space;
    VOLMIN=VMIN;
    FILLFACTOR=FFACTOR;
    if(rnode == NULL)
      {
	try
	  {
	    rnode=new OcTreeNode;
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
    bool smallVOL=VOL <= VOLMIN;
    double ff=(double)(rnode->ppoints.size())/(double)(VOL);
    bool veryFULL= ff >= FILLFACTOR;
    if(smallVOL || veryFULL)
      {
	Point* pFAKE=0;
	for(int ni=rnode->ppoints.size();ni<VOL;ni++)
	  rnode->ppoints.push_back(pFAKE);
      }
    rnode->full=rnode->ppoints.size() == VOL;
    if(rnode->full)
      {
	fullnodes++;
	// cerr << " LOAD NODES0 " << RANK << " " << nnodes << " " << fullnodes << " " << vol << " " << pPOINTS.size() << endl;
	return;
      }
    for(int ni=0;ni<2;ni++)
      if(rnode->kids[ni] == NULL)
	{
	  try
	    {
	      rnode->kids[ni]= new OcTreeNode;
	    }
	  catch(bad_alloc& ba)
	    {
	      cerr << " badd root kid " << ba.what() << " " << ni << " " << RANK << endl;
	      exit(0);
	    }
	}
    for(int corner=0;corner<2;corner++)
      {
	rnode->kids[corner]->dir=2;
	LoadOcTree(corner,rnode);
      }
    // cerr << " LOADR " << RANK << " " << nnodes << " " << fullnodes << " " << pPOINTS.size() << " " << vol << " ";
    // cerr << BOX[0] << " " << BOX[1] << " " << BOX[2] << " " << BOX[3] << " " << BOX[4] << " " << BOX[5] << endl;
    nnodes=0;
    fullnodes=0;
  }
  void OcTree::LoadOcTree(int cornera,OcTreeNode* pnode)
  {
    nnodes++;
    OcTreeNode* knode=pnode->kids[cornera];
    knode->kids.assign(2,NULL);
    knode->box=pnode->box;
    int dir=knode->dir;
    int DIR2=2*dir;
    knode->box[cornera % 2 == 0 ? DIR2+1:DIR2]=(pnode->box[DIR2]+pnode->box[DIR2+1])/2;
    // cerr << " LOADP " << RANK << " " << vol << " " << DIR;
    // cerr << " " << knode->box[0] << " " << knode->box[1] << " " << knode->box[2] << " " << knode->box[3] << " " << knode->box[4] << " " << knode->box[5] << endl;
    knode->ppoints.clear();
    vector <int>pos(3);
    int MINX=Misc::pow(2,29);
    int MAXX=-MINX;
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
	    MINX=min(MINX,pos[dir]);
	    MAXX=max(MAXX,pos[dir]);
	  }
      }
    knode->box[DIR2]=MINX;
    knode->box[DIR2+1]=MAXX+1;
    int vol=(knode->box[1]-knode->box[0]);
    vol*=(knode->box[3]-knode->box[2]);
    vol*=(knode->box[5]-knode->box[4]);
    knode->empty=vol == 0 || knode->ppoints.empty();
    if(knode->empty)
      {
	knode->full=false;
	return;
      }
    bool smallVOL=vol <= VOLMIN;
    double ff=(double)(knode->ppoints.size())/(double)(vol);
    bool veryFULL= ff >= FILLFACTOR;
    if(smallVOL || veryFULL)
      {
	Point* pFAKE=0;
	for(int ni=knode->ppoints.size();ni<vol;ni++)
	  knode->ppoints.push_back(pFAKE);
      }
    knode->full=vol != 0 && knode->ppoints.size() == vol;
    if(knode->full)
      {
	fullnodes++;
      }
    if(knode->full)
      return;
    for(int ni=0;ni<2;ni++)
      if(knode->kids[ni] == NULL)
	{
	  try
	    {
	      knode->kids[ni]= new OcTreeNode;
	    }
	  catch(bad_alloc& ba)
	    {
	      cerr << " badd node kid " << ba.what() << " " << ni << " " << RANK << endl;
	      exit(0);
	    }
	}
    for(int cornerb=0;cornerb<2;cornerb++)
      {
	knode->kids[cornerb]->dir=(knode->dir+1)%3;
	LoadOcTree(cornerb,knode);
      }
  }   
  void OcTree::DestroyOcTree()
  {
    DestroyOcTree(rnode);
  }
  void OcTree::DestroyOcTree(OcTreeNode* pnode)
  {
    if(pnode != NULL)
      {
	for(int k=0;k<2;k++)
	  DestroyOcTree(pnode->kids[k]);
	delete pnode;
      }
  }
  void OcTree::CollectBoxes(vector < vector<int> >& SBoxes)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	// int nb=SBoxes.size();
	// SBoxes.resize(nb+1);
	// SBoxes[nb].assign(rnode->box.begin(),rnode->box.end());
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=rnode->box;
	return;
      }
    for(int k=0;k<2;k++)
      CollectBoxes(SBoxes,rnode->kids[k]);
  }
  void OcTree::CollectBoxes(vector< vector<int> >& SBoxes,OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	// int nb=SBoxes.size();
	// SBoxes.resize(nb+1);
	// SBoxes[nb].assign(pnode->box.begin(),pnode->box.end());
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=pnode->box;
	return;
      }
    for(int k=0;k<2;k++)
      CollectBoxes(SBoxes,pnode->kids[k]);
  }
  void OcTree::CollectPoints(vector < vector<Point*> >& SPoints)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	// int np=SPoints.size();
	// SPoints.resize(np+1);
	// SPoints[np].assign(rnode->ppoints.begin(),rnode->ppoints.end());
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    for(int k=0;k<2;k++)
      CollectPoints(SPoints,rnode->kids[k]);
  }
  void OcTree::CollectPoints(vector< vector<Point*> >& SPoints,OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	// int np=SPoints.size();
	// SPoints.resize(np+1);
	// SPoints[np].assign(pnode->ppoints.begin(),pnode->ppoints.end());
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    for(int k=0;k<2;k++)
      CollectPoints(SPoints,pnode->kids[k]);
  }
  void OcTree::CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	// int nb=SBoxes.size();
	// SBoxes.resize(nb+1);
	// SBoxes[nb].assign(rnode->box.begin(),rnode->box.end());
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=rnode->box;
	// int np=SPoints.size();
	// SPoints.resize(np+1);
	// SPoints[np].assign(rnode->ppoints.begin(),rnode->ppoints.end());
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    for(int k=0;k<2;k++)
      CollectBoxesPoints(SBoxes,SPoints,rnode->kids[k]);
  }
  void OcTree::CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints,OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	// int nb=SBoxes.size();
	// SBoxes.resize(nb+1);
	// SBoxes[nb].assign(pnode->box.begin(),pnode->box.end());
	SBoxes.resize(SBoxes.size()+1);
	SBoxes.back()=pnode->box;
	// int np=SPoints.size();
	// SPoints.resize(np+1);
	// SPoints[np].assign(pnode->ppoints.begin(),pnode->ppoints.end());
	// np++;
	SPoints.resize(SPoints.size()+1);
	SPoints.back().assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    for(int k=0;k<2;k++)
      CollectBoxesPoints(SBoxes,SPoints,pnode->kids[k]);
  }
  void OcTree::DisplayTree(int& TOT,int& NB)
  {
    if(rnode != NULL)
      {
	int vol=(rnode->box[1]-rnode->box[0]);
	vol*=(rnode->box[3]-rnode->box[2]);
	vol*=(rnode->box[5]-rnode->box[4]);
	if(rnode->full)
	  {
	    assert(rnode->ppoints.size() == vol);
	    TOT+=rnode->ppoints.size();
	    NB++;
	    return;
	  }
	else
	  assert(rnode->ppoints.size() < vol);
	DisplayTree(rnode,TOT,NB);
      }
  }
  void OcTree::DisplayTree(OcTreeNode* pnode,int& TOT,int& NB)
  {
    if(pnode != NULL)
      {
	int vol=(pnode->box[1]-pnode->box[0]);
	vol*=(pnode->box[3]-pnode->box[2]);
	vol*=(pnode->box[5]-pnode->box[4]);
	if(pnode->full)
	  {
	    assert(pnode->ppoints.size() == vol);
	    TOT+=pnode->ppoints.size();
	    NB++;
	    return;
	  }
	else
	  if(vol > 0)
	    assert(pnode->ppoints.size() < vol);
	for(int k=0;k<2;k++)
	  DisplayTree(pnode->kids[k],TOT,NB);
      }
  }
  void OcTree::Traverse()
  {
    if(rnode != NULL)
      Traverse(rnode);
  }
  void OcTree::Traverse(OcTreeNode* pnode)
  {
    if(pnode != NULL)
      for(int k=0;k<2;k++)
	Traverse(pnode->kids[k]);
  }
}
