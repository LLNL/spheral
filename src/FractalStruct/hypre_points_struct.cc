#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_struct(Fractal_Memory& mem,vector <Group*>& groups,
			   vector < vector <Point*> >& hypre_points,bool buffer_groups,int level)
  {
    int RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    bool Ranky = RANK == 21;
    Ranky=true;
    vector <int>pos(3);
    vector <int> BOX=mem.BoxesLev[mem.p_mess->FractalRank][level];
    hypre_points.clear();
    int ng=0;
//     for(vector <Group*>::const_iterator group_itr=groups.begin();group_itr!=groups.end();group_itr++)
    for(Group* pgroup : groups)
      {
// 	Group* pgroup=*group_itr;
 	if(buffer_groups == pgroup->get_buffer_group())
// 	if(true)
	  {
	    hypre_points.resize(ng+1);
// 	    for(vector<Point*>::const_iterator point_itr=pgroup->list_points.begin();point_itr !=pgroup->list_points.end();++point_itr)
	    for(Point* p : pgroup->list_points)
	      {
// 		Point* p=*point_itr;
		p->get_pos_point(pos);
		if(p->get_inside() && vector_in_box(pos,BOX))
		  hypre_points[ng].push_back(p);
	      }
	    if(!hypre_points[ng].empty())
	      sort3_list(hypre_points[ng],0);
	    if(Ranky)
	      {
		cerr << " GROUPA " << pgroup->list_points.size() << endl;
		cerr << " GROUPS " << ng << " " << (hypre_points[ng]).size() << endl;
	      }
	    ng++;
	  }
      }
    cerr << " GroupTotal " << RANK << " " << ng << endl;
  }
}
