#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void sort_3(Fractal& fractal,Group& group)
  {
    Point* haha=0;
    vector <Point*> which(27,haha);
    vector <Point*> ud0(3,haha);
    int cx,cy,cz;
    for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point*p=*point_itr;
	which[p->get_real_pointer()]=p;
	p->set_point_ud(ud0);
      }
    for(int counter=0;counter < 26;counter++)
      {
	cx=counter % 3;
	cy=counter/3 % 3;
	cz=counter/9;
	if(cx < 2)
	  which[counter]->set_point_up_x(which[counter+1]);
	if(cy < 2)
	  which[counter]->set_point_up_y(which[counter+3]);
	if(cz < 2)
	  which[counter]->set_point_up_z(which[counter+9]);
      }
  }
  /*
   void sort_3(Fractal& fractal,Group& group,const int& what1,const bool& re_order,const bool& remove_duplicates)
   {
     const int length=fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max());
     const int spacing=Misc::pow(2,fractal.get_level_max()-group.get_level());
     vector <list <Point*> > list_a;
     sort_list_1(group.list_points,list_a,length,spacing,what1,1);
       cout << " start of 3, list size " << list_a.size() << endl;
     int maxlists=fractal.get_grid_length()*Misc::pow(2,group.get_level());
     int total_c=-1;
     Point* haha=0;
     vector <Point*> list_c(maxlists,haha);
     bool period=fractal.get_periodic();
     int what1t2m1=2*what1-1;
     if(re_order)
       group.list_points.clear();
     int total_a=list_a.size();
     for(int t_a=0;t_a < total_a;++t_a)
       {
         if(!list_a[t_a].empty())
   	{
   	  vector <list <Point*> > list_b;
   	  sort_list_1(list_a[t_a],list_b,length,spacing,what1,2);
   	  int total_b=list_b.size();
   	  for(int t_b=0;t_b < total_b;++t_b)
   	    {
   	      if(!list_b[t_b].empty())
   		{
   		  sort_list_1(list_b[t_b],list_c,remove_duplicates,total_c,spacing,what1,3); 
   		  if(total_c > 1) 
   		    {
   		      for(int t_c=0;t_c < total_c-1;++t_c)
   			{
   			  if(list_c[t_c] != 0 && list_c[t_c+1] != 0)
   			    list_c[t_c]->set_point_ud(list_c[t_c+1],what1t2m1);
   			}
   		      if(period && total_c == maxlists)
   			list_c[total_c-1]->set_point_ud(list_c[0],what1t2m1);
   		    }
   		  if(re_order && total_c >0)
   		    {
   		      for(int t_c=0;t_c < total_c;++t_c)
   			{
   			  if(list_c[t_c] != 0)
   			    group.list_points.push_back(list_c[t_c]);
   			}
   		    }
   		}
   	    }
   	}
       }
   }
   void sort_list_1(list <Point*>& list_points,vector <list <Point*> >& list_a,
   		  const int& length,const int& spacing,const int& what1,const int& what2)
   {
     if(list_points.empty()) return;
     int maxx=-INT_MAX;
     int minn=INT_MAX;
     int what=2-(what2-what1+3)%3;
     int total=-1;
     for(list <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
       {
         int poss=(*point_itr)->get_pos_point(what);
         maxx=max(maxx,poss);
         minn=min(minn,poss);
       }
     total=(maxx-minn)/spacing+1;
     list_a.resize(total);
     for(list <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
       {
         int t=((*point_itr)->get_pos_point(what)-minn)/spacing;
         assert(t >=0 && t < total);
         list_a[t].push_back(*point_itr);
       }
   }
   void sort_list_1(list <Point*>& list_points,vector <Point*>& list_a,const bool& remove,
   		  int& total,const int& spacing,const int& what1,const int& what2)
   {
       cout << "correct sort" << endl;
     if(list_points.empty()) return;
     int maxx=-INT_MAX;
     int minn=INT_MAX;
     int what=2-(what2-what1+3)%3;
     for(list <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
       {
         int poss=(*point_itr)->get_pos_point(what);
         maxx=max(maxx,poss);
         minn=min(minn,poss);
       }
     total=(maxx-minn)/spacing+1;
     Point* pp=0;
       list_a.assign(total,pp);
     Misc::my_assign(list_a,0,total,pp);
     for(list <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
       {
         int t=((*point_itr)->get_pos_point(what)-minn)/spacing;
         assert(t >=0 && t < total);
         if(remove && list_a[t] !=0)
   	{
   	  if(list_a[t]->get_real_pointer() > (*point_itr)->get_real_pointer())
   	    list_a[t]=*point_itr;
   	}
         else
   	list_a[t]=*point_itr;
       }
   }
*/
}
