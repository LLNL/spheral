#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::pot_solver="sor";
  //
  void sor(Group& group, Fractal& fractal,vector <Point*>& list_edge,const int& dir)
  {
    ofstream& FileSor=fractal.p_file->DUMPS;
    double rj=group.get_rjac();
    double g_c=group.get_force_const();
    //
    double anormf=0.0;
    int n_points=0;
    //
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(point.get_inside())
	  {
	    n_points++;
	    anormf+=abs(point.get_potential_point());
	  }
      }
    double anorm2f=3.0*anormf/(double)n_points;
    int maxits_real=fractal.get_maxits();
    if(n_points==1)maxits_real=1;
    double omega=1.0;
    int n=0;
    double anorm=0.0;
    double anorm2=0.0;
    for(n=0;n< maxits_real;++n)
      {
	anorm=0.0;
	anorm2=-1.0;
	for(int ipass=0;ipass < 2;++ipass)
	  {
	    double omega_6=omega/6.0;
	    bool ipass1=ipass==1;
	    for (vector<Point*>::const_iterator  point_itr=list_edge.begin();point_itr !=list_edge.end();++point_itr)
	      {
		Point* p_point=*point_itr;
		if(ipass1)p_point=p_point->get_point_ud_0(dir,2);
		bool go_on=it_is_inside(p_point); 
		while(go_on)
		  {
		    double resid=p_point->laplacian()-g_c*p_point->get_density_point();
		    anorm+=abs(resid);
		    anorm2=max(anorm2,abs(resid));
		    p_point->add_potential_point(omega_6*resid);
		    p_point=p_point->get_point_ud_0(dir,3);
		    p_point=p_point->get_point_ud(dir);
		    go_on=it_is_inside(p_point);
		  }
	      }
	    if(n==0 && ipass ==0)
	      omega=1.0/(1.0-0.5*rj*rj);
	    else
	      omega=1.0/(1.0-0.25*rj*rj*omega);
	  }
	if(anorm < fractal.get_epsilon_sor()*anormf && anorm2 < fractal.get_epsilon_sor()*anorm2f) 
	  {
// 	    FileSor << n << "\t " 
// 		    << n_points << "\t " 
// 		    << scientific
// 		    << anorm << "\t " 
// 		    << anormf << "  " 
// 		    << anorm2 << "\t " 
// 		    << anorm2f  << "\t " 
// 		    << rj << "\t " 
// 		    << fractal.get_epsilon_sor() << "\t"
// 		    << group.get_level() 
// 		    << "\n";
	    return;
	  }
      }
    if(maxits_real ==1) return;
    FileSor << " not converged " << "\n";
    /*
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	point.dumpp();
      }
    */
    FileSor << n << "\t " 
	    << n_points << "\t " 
	    << scientific
	    << anorm << "\t " 
	    << anormf << "  " 
	    << anorm2 << "\t " 
	    << anorm2f  << "\t " 
	    << rj << "\t " 
	    << fractal.get_epsilon_sor() << "\t"
	    << group.get_level() 
	    << "\n";
    FileSor << " not converged " << "\n";
  }
  bool it_is_inside(Point* p_point)
  {
    return p_point != 0 && p_point->get_inside();
  }
}
