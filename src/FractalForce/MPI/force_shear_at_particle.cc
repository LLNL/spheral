#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_shear_at_particle(Fractal_Memory& fractal_memory,Fractal& fractal,Misc& misc)
  { 
    vector <double> dens(8);
    vector <double> weights(8);
    vector <double> f_xx(8);
    vector <double> f_xy(8);
    vector <double> f_xz(8);
    vector <double> f_yy(8);
    vector <double> f_yz(8);
    vector <double> f_zz(8);
    vector <double> shear(6);
    vector <double> shear2(6);
    shear2.assign(6,0.0);
    const  double scale=(double)(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
    double fourpiinv=1.0/(16.0*atan(1.0));
    //
    for(int level=0;level <= fractal.get_level_max();level++)
      {
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group* p_group=*group_itr;
	    Group& group=*p_group;
	    double d_inv=pow(2.0,group.get_level()-fractal.get_level_max());
	    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point& point=**point_itr;
		//		if(point.get_buffer_point()) continue;
		bool not_yet=true;
		//
		//	cout << "doing point " << *point_itr << endl;
		for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
		  {
		    Particle& particle=**particle_itr;
		    // cout << "doing particle " << *particle_itr << endl;
		    // cout << "p group " << p_group << " " << particle.get_p_highest_level_group() << endl;
		    if(p_group == particle.get_p_highest_level_group())
		      {
			if(not_yet)
			  {
			    // cout << " not yet a " << &point << endl;
			    point.get_field_shear_values(f_xx,f_xy,f_xz,f_yy,f_yz,f_zz);
			    not_yet=false;
			    // cout << " not yet b " << endl;
			  }
			//
			double d_x,d_y,d_z;
			vector <double> pos(3);
			// cout << "pos a" << endl;
			particle.get_pos(pos);
			// cout << "pos b" << endl;
			point.get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
			// cout << "delta a" << endl;
			Misc::set_weights(weights,d_x,d_y,d_z);
			Misc::sum_prod<double>(0,7,1,shear,weights,f_xx,f_xy,f_xz,f_yy,f_yz,f_zz);
			if(fractal_memory.calc_density_particle)
			  {
			    // cout << "shear densa " << endl;
			    particle.set_density(-(shear[0]+shear[3]+shear[5])*fourpiinv);
			    // cout << "shear densb " << endl;
			  }
			if(fractal_memory.calc_shear)
			  {
			    for(int i=0;i<6;i++)
			      shear2[i]+=shear[i]*shear[i];
			    double r=0.0;
			    max_predict(fractal_memory,fractal,shear,r);
			    if(fractal_memory.start_up) 
			      {
				// cout << "sheara " << r << endl;
				particle.set_rad_max(r);
				// cout << "shearb " << r << endl;
			      }
			  }
		      }
		  }
	      }
	  }
      }

    for(int level=0;level <= fractal.get_level_max();level++)
      {
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      (*point_itr)->force_shear_point_zero();
	  }
      }


    if(fractal_memory.calc_shear)
      {
	double spam=(double)fractal.get_number_particles();
	cout << " shear2 " ;
	cout << shear2[0]/spam << " " ;
	cout << shear2[1]/spam << " " ;
	cout << shear2[2]/spam << " " ;
	cout << shear2[3]/spam << " " ;
	cout << shear2[4]/spam << " " ;
	cout << shear2[5]/spam << " " ;
	cout << endl;
	//    assert(0);
      }
  }
}
