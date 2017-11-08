#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void clean_groups(Fractal_Memory& mem)
  {
    int totalsizeA=0;
    int totalsizeB=0;
    int totalcapA=0;
    int totalcapB=0;
    for(int level=0;level<=mem.p_fractal->get_level_max();level++)
      {
	for(auto &pg : mem.all_groups[level])
	  {
	    {
	      totalsizeA+=pg->p_list_really_high.size();
	      totalcapA+=pg->p_list_really_high.capacity();
	      vector<Point*>holy_grail;
	      pg->p_list_really_high.swap(holy_grail);
	      totalsizeB+=pg->p_list_really_high.size();
	      totalcapB+=pg->p_list_really_high.capacity();
	    }
	    {
	      totalsizeA+=pg->list_high_points.size();
	      totalcapA+=pg->list_high_points.capacity();
	      vector<Point*>Icontinentia_Maxima;
	      pg->list_high_points.swap(Icontinentia_Maxima);	    
	    }
	    {
	      totalsizeA+=pg->list_high.size();
	      totalcapA+=pg->list_high.capacity();
	      vector<Point*>Biggus_Dickus;
	      pg->list_high.swap(Biggus_Dickus);
	      totalsizeB+=pg->list_high.size();
	      totalcapB+=pg->list_high.capacity();
	    }
	    {
	      totalsizeA+=pg->list_high_groups.size();
	      totalcapA+=pg->list_high_groups.capacity();
	      vector<Group*>Castle_Anthrax;
	      pg->list_high_groups.swap(Castle_Anthrax);
	      totalsizeB+=pg->list_high_groups.size();
	      totalcapB+=pg->list_high_groups.capacity();
	    }
	    {
	      totalsizeA+=pg->head_number.size();
	      totalcapA+=pg->head_number.capacity();
	      vector<int>Sir_Robin;
	      pg->head_number.swap(Sir_Robin);
	      totalsizeB+=pg->head_number.size();
	      totalcapB+=pg->head_number.capacity();
	    }
	    {
	      totalsizeA+=pg->list_pair_1.size();
	      totalcapA+=pg->list_pair_1.capacity();
	      vector<int>Judith;
	      pg->list_pair_1.swap(Judith);
	      totalsizeB+=pg->list_pair_1.size();
	      totalcapB+=pg->list_pair_1.capacity();
	    }
	    {
	      totalsizeA+=pg->list_pair_2.size();
	      totalcapA+=pg->list_pair_2.capacity();
	      vector<int>black_knight;
	      pg->list_pair_2.swap(black_knight);
	      totalsizeB+=pg->list_pair_2.size();
	      totalcapB+=pg->list_pair_2.capacity();
	    }
	  }
      }
    ofstream& FHT=mem.p_file->DUMPS;
    FHT << "Clean Groups " << mem.steps << " " << totalsizeA << " " << totalsizeB << " ";
    FHT << totalcapA << " " << totalcapB << "\n";
  }
}
