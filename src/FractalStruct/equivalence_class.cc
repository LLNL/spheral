#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void equivalence_class(Group& group)
  {
    group.head_number.resize(group.get_number_high_points());
    for(int k=0; k < group.get_number_high_points(); ++k)
      group.head_number[k]=k;
    for(int i=0; i < group.get_number_high_pairs(); ++i)
      {
	int j=group.list_pair_1[i];
	while(group.head_number[j]  != j)
	  j=group.head_number[j];
	int k=group.list_pair_2[i];
	while(group.head_number[k] != k)
	  k=group.head_number[k];
	if(j != k) 
	  group.head_number[j]=k;
      }
    for(int j=0;j < group.get_number_high_points(); ++j)
      while(group.head_number[j] != group.head_number[group.head_number[j]])
	group.head_number[j]=group.head_number[group.head_number[j]];
    group.list_pair_1.clear();
    group.list_pair_2.clear();
  }
}
