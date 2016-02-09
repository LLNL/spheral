#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  OcTreeNode::OcTreeNode()
  {
    full=false;
    box=NULL;
    kids=NULL;
  }
  OcTreeNode::~OcTreeNode()
  {
    if(box != NULL)
      delete [] box;
    if(kids != NULL)
      delete [] kids;
  }
  OcTree::OcTree()
  {
    rnode=NULL;
    VolumeFactor=1;
  }
  OcTree::~OcTree()
  {
    if(rnode != NULL)
      delete rnode;
  }
  void OcTree::LoadOcTree(int* BOX,vector <Point*>& pPOINTS,int spacing=1)
  {
    VolumeFactor=spacing*spacing*spacing;
    if(rnode == NULL)
      rnode=new OcTreeNode;
    if(rnode->box == NULL)
      rnode->box=new int[6];
    if(rnode->kids == NULL)
      rnode->kids= new OcTreeNode*[8];
    std::copy(BOX,BOX+6,rnode->box);
    rnode->ppoints.assign(pPOINTS.begin(),pPOINTS.end());
    int vol=(BOX[1]-BOX[0])*(BOX[3]-BOX[2])*(BOX[5]-BOX[4]);
    rnode->full=rnode->ppoints.size()*VolumeFactor == vol;
    if(rnode->ppoints.empty() || rnode->full);
      return;
    for(int corner=0;corner<8;corner++)
      LoadOcTree(corner,rnode);
  }
  void OcTree::LoadOcTree(int cornera,OcTreeNode* pnode)
  {
    OcTreeNode* knode=pnode->kids[cornera];
    if(knode->box == NULL)
      knode->box=new int[6];
    std::copy(pnode->box,pnode->box+6,knode->box);
    knode->box[cornera % 2 == 0 ? 1:0]=(pnode->box[0]+pnode->box[1])/2;
    knode->box[(cornera/2) % 2 == 0 ? 3:2]=(pnode->box[2]+pnode->box[3])/2;
    knode->box[cornera/4 == 0 ? 5:4]=(pnode->box[4]+pnode->box[5])/2;
    int vol=(knode->box[1]-knode->box[0]);
    vol*=(knode->box[3]-knode->box[2]);
    vol*=(knode->box[5]-knode->box[4]);
    knode->ppoints.clear();
    vector <int>pos(3);
    for(list<Point*>::iterator itp=pnode->ppoints.begin();itp!=pnode->ppoints.end();itp++)
      {
	(*itp)->get_pos_point(pos);
	if(pos[0] < knode->box[0] || pos[1] < knode->box[2] || pos[2] < knode->box[4] ||
	   pos[0] >= knode->box[1] || pos[1] >= knode->box[3] || pos[2] >= knode->box[5])
	  continue;
	knode->ppoints.push_back(*itp);
	itp=pnode->ppoints.erase(itp);
	if(itp==pnode->ppoints.end())
	  break;
      }
    knode->full=vol != 0 && knode->ppoints.size()*VolumeFactor == vol;
    if(knode->ppoints.empty() || knode->full)
      return;
    if(knode->kids == NULL)
      knode->kids=new OcTreeNode*[8];
    for(int cornerb=0;cornerb<8;cornerb++)
      LoadOcTree(cornerb,knode);
  }   
  void OcTree::DestroyOcTree()
  {
    DestroyOcTree(rnode);
    delete rnode;
  }
  void OcTree::DestroyOcTree(OcTreeNode* pnode)
  {
    if(pnode == NULL)
      return;
    for(int k=0;k<8;k++)
      DestroyOcTree(pnode->kids[k]);
    delete pnode;
  }
  void OcTree::CollectBoxes(vector < vector<int> >& SSBoxes)
  {
    if(rnode->ppoints.empty())
      return;
    if(rnode->full)
      {
	SSBoxes.resize(1);
	SSBoxes[0].resize(6);
	std::copy(rnode->box,rnode->box+6,SSBoxes[0].begin());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxes(SSBoxes,rnode->kids[k]);
  }
  void OcTree::CollectBoxes(vector< vector<int> >& SSBoxes,OcTreeNode* pnode)
  {
    if(pnode->ppoints.empty())
      return;
    if(pnode->full)
      {
	SSBoxes.resize(SSBoxes.size()+1);
	SSBoxes.back().resize(6);
	std::copy(pnode->box,pnode->box+6,SSBoxes.back().begin());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxes(SSBoxes,pnode->kids[k]);
  }
  void OcTree::CollectPoints(vector < vector<Point*> >& SSPoints)
  {
    if(rnode->ppoints.empty())
      return;
    if(rnode->full)
      {
	SSPoints.resize(1);
	SSPoints[0].assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectPoints(SSPoints,rnode->kids[k]);
  }
  void OcTree::CollectPoints(vector< vector<Point*> >& SSPoints,OcTreeNode* pnode)
  {
    if(pnode->ppoints.empty())
      return;
    if(pnode->full)
      {
	SSPoints.resize(SSPoints.size()+1);
	(SSPoints.back()).assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectPoints(SSPoints,pnode->kids[k]);
  }
  void OcTree::CollectBoxesPoints(vector < vector<int> >& SSBoxes,vector < vector<Point*> >& SSPoints)
  {
    if(rnode->ppoints.empty())
      return;
    if(rnode->full)
      {
	SSBoxes.resize(1);
	std::copy(rnode->box,rnode->box+6,SSBoxes[0].begin());
	SSPoints.resize(1);
	SSPoints[0].assign(rnode->ppoints.begin(),rnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxesPoints(SSBoxes,SSPoints,rnode->kids[k]);
  }
  void OcTree::CollectBoxesPoints(vector < vector<int> >& SSBoxes,vector < vector<Point*> >& SSPoints,OcTreeNode* pnode)
  {
    if(pnode->ppoints.empty())
      return;
    if(pnode->full)
      {
	SSBoxes.resize(SSBoxes.size()+1);
	std::copy(pnode->box,pnode->box+6,(SSBoxes.back()).begin());
	SSPoints.resize(SSPoints.size()+1);
	SSPoints.back().assign(pnode->ppoints.begin(),pnode->ppoints.end());
	return;
      }
    for(int k=0;k<8;k++)
      CollectBoxesPoints(SSBoxes,SSPoints,pnode->kids[k]);
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
}
