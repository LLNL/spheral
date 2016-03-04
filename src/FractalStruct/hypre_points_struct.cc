#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_struct(Fractal_Memory& mem,vector <Group*>& groups,
			   vector < vector <Point*> >& hypre_points,bool buffer_groups,int level)
  {
//     int RANK=-1;
//     MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
//     bool Ranky = RANK == 21;
//     Ranky=true;
    vector <int>pos(3);
    vector <int> BOX=mem.BoxesLev[mem.p_mess->FractalRank][level];
    hypre_points.clear();
    int ng=0;
    for(Group* &pgroup : groups)
      {
 	if(buffer_groups == pgroup->get_buffer_group())
	  {
	    hypre_points.resize(ng+1);
	    for(Point* &p : pgroup->list_points)
	      {
		p->get_pos_point(pos);
		if(p->get_inside() && vector_in_box(pos,BOX))
		  hypre_points[ng].push_back(p);
	      }
	    if(!hypre_points[ng].empty())
	      sort3_list(hypre_points[ng],0);
// 	    if(Ranky)
// 	      {
// 		cerr << " GROUPA " << RANK << " " << pgroup->list_points.size();
// 		cerr << " GROUPS " << ng << " " << hypre_points[ng].size() << endl;
// 	      }
	    ng++;
	  }
      }
    cerr << " GroupTotal " << RANK << " " << ng << endl;
  }
}
