#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  OcTreeNode::OcTreeNode():
    full(false)
  {}
  OcTreeNode::~OcTreeNode()
  {}
  OcTree::OcTree():
    nnodes(0),
    fullnodes(0),
    rnode(NULL),
    spacing(1)
  {
    RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
  }
  OcTree::~OcTree()
  {
    DestroyOcTree();
  }
  void OcTree::LoadOcTree(vector<int>& BOX,vector <Point*>& pPOINTS,int space=1):
    spacing(space)
  {
    nnodes++;
    // spacing=space;
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
    rnode->kids.assign(8,NULL);
    rnode->box=BOX;
    rnode->ppoints.assign(pPOINTS.begin(),pPOINTS.end());
    int vol=(BOX[1]-BOX[0])*(BOX[3]-BOX[2])*(BOX[5]-BOX[4]);
    rnode->full=rnode->ppoints.size() == vol;
    if(rnode->full)
      {
	fullnodes++;
	// cerr << " LOAD NODES0 " << RANK << " " << nnodes << " " << fullnodes << " " << vol << " " << pPOINTS.size() << endl;
      }
    if(rnode->ppoints.empty() || rnode->full)
      return;
    for(int ni=0;ni<8;ni++)
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
    for(int corner=0;corner<8;corner++)
      LoadOcTree(corner,rnode);
    // cerr << " LOAD NODES0 " << RANK << " " << nnodes << " " << fullnodes << " " << pPOINTS.size() << " " << vol << " ";
    // cerr << BOX[0] << " " << BOX[1] << " " << BOX[2] << " " << BOX[3] << " " << BOX[4] << " " << BOX[5] << endl;
    nnodes=0;
    fullnodes=0;
  }
  void OcTree::LoadOcTree(int cornera,OcTreeNode* pnode)
  {
    nnodes++;
    OcTreeNode* knode=pnode->kids[cornera];
    knode->kids.assign(8,NULL);
    knode->box=pnode->box;
    knode->box[cornera % 2 == 0 ? 1:0]=(pnode->box[0]+pnode->box[1])/2;
    knode->box[(cornera/2) % 2 == 0 ? 3:2]=(pnode->box[2]+pnode->box[3])/2;
    knode->box[cornera/4 == 0 ? 5:4]=(pnode->box[4]+pnode->box[5])/2;
    int vol=(knode->box[1]-knode->box[0]);
    vol*=(knode->box[3]-knode->box[2]);
    vol*=(knode->box[5]-knode->box[4]);
    knode->ppoints.clear();
    vector <int>pos(3);
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
	  }
      }
    knode->full=vol != 0 && knode->ppoints.size() == vol;
    if(knode->full)
      {
	fullnodes++;
      }
    if(knode->ppoints.empty() || knode->full)
      return;
    for(int ni=0;ni<8;ni++)
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
    for(int cornerb=0;cornerb<8;cornerb++)
      LoadOcTree(cornerb,knode);
  }   
  void OcTree::DestroyOcTree()
  {
    DestroyOcTree(rnode);
  }
  void OcTree::DestroyOcTree(OcTreeNode* pnode)
  {
    if(pnode != NULL)
      {
	for(int k=0;k<8;k++)
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
	int nb=SBoxes.size();
	SBoxes.resize(nb+1);
	SBoxes[nb].assign(rnode->box.begin(),rnode->box.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxes(SBoxes,rnode->kids[k]);
  }
  void OcTree::CollectBoxes(vector< vector<int> >& SBoxes,OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	int nb=SBoxes.size();
	SBoxes.resize(nb+1);
	SBoxes[nb].assign(pnode->box.begin(),pnode->box.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxes(SBoxes,pnode->kids[k]);
  }
  void OcTree::CollectPoints(vector < vector<Point*> >& SPoints)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	int np=SPoints.size();
	SPoints.resize(np+1);
	SPoints[np].assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectPoints(SPoints,rnode->kids[k]);
  }
  void OcTree::CollectPoints(vector< vector<Point*> >& SPoints,OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	int np=SPoints.size();
	SPoints.resize(np+1);
	SPoints[np].assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectPoints(SPoints,pnode->kids[k]);
  }
  void OcTree::CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    if(rnode == NULL)
      return;
    if(rnode->full)
      {
	int nb=SBoxes.size();
	SBoxes.resize(nb+1);
	SBoxes[nb].assign(rnode->box.begin(),rnode->box.end());
	int np=SPoints.size();
	SPoints.resize(np+1);
	SPoints[np].assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxesPoints(SBoxes,SPoints,rnode->kids[k]);
  }
  void OcTree::CollectBoxesPoints(vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints,OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    if(pnode->full)
      {
	int nb=SBoxes.size();
	SBoxes.resize(nb+1);
	// 	SBoxes.back()=pnode->box;
	SBoxes[nb].assign(pnode->box.begin(),pnode->box.end());
	nb++;
	int np=SPoints.size();
	SPoints.resize(np+1);
	SPoints[np].assign(pnode->ppoints.begin(),pnode->ppoints.end());
	np++;
	return;
      }
    for(int k=0;k<8;k++)
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
	for(int k=0;k<8;k++)
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
      for(int k=0;k<8;k++)
	Traverse(pnode->kids[k]);
  }
  // void OcTree::Consolidate()
  // {
  //   if(rnode != NULL && !rnode->full && !rnode->ppoints.empty())
  //     Consolidate(rnode);
  // }
  // void OcTree::Consolidate(OcTreeNode* pnode)
  // {
  //   if(pnode == NULL)
  //     return;
  //   vector <bool>fully(8,false);
  //   for(int k=0;k<8;k++)
  //     if(pnode->kids[k] != NULL)
  // 	fully[k]=pnode->kids[k]->full;
  //   if(fully[0] && fully[1] && fully[2] &&fully[3])
  //     {
  // 	for(int b : {1,3,5})
  // 	  pnode->kids[0]->box[b]=pnode->kids[3]->box[b];
  // 	for(int k : {1,2,3})
  // 	  {
  // 	    pnode->kids[k]->full=false;
  // 	    std::copy(pnode->kids[k]->ppoints.begin(),pnode->kids[k]->ppoints.begin(),pnode->kids[0]->ppoints.end());
  // 	    pnode->kids[k]->ppoints.clear();
  // 	  }
  //     }
  //   for(int k=0;k<8;k++)
  //     Consolidate(pnode->kids[k]);
  // }
}
