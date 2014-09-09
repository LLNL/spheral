#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <class M, class F>  int split_particle(M& mem,F& frac,const double& x0,const double& y0,const double& z0,
						  int& count,const double& m,const int& split_to,const bool& gen_part)
  {
    Particle* par=0;
    int iphase=3;
    if(mem.do_vel) iphase=6;
    int jfield=6;
    //    if(mem.calc_density_particle) jfield=5;
    //    if(mem.calc_shear) jfield=12;
    int many_split=1;
    int level_mask=0;
    many_split=max(split_to,1);
    for(int i=0;i<mem.splits;i++)
      {
	bool inside=false;
	if(mem.splits_many[i] > many_split)
	  {
	    double dx=x0-mem.splits_center_x[i];
	    double dy=y0-mem.splits_center_y[i];
	    double dz=z0-mem.splits_center_z[i];
	    inside=false;
	    if(mem.periodic)
	      {
		dx=min(min(abs(dx),abs(dx+1.0)),abs(dx-1.0));
		dy=min(min(abs(dy),abs(dy+1.0)),abs(dy-1.0));
		dz=min(min(abs(dz),abs(dz+1.0)),abs(dz-1.0));
	      }
	    if(mem.splits_square[i])
	      {
		inside=
		  abs(dx/mem.splits_rad_x[i]) <= 1.0 &&
		  abs(dy/mem.splits_rad_y[i]) <= 1.0 &&
		  abs(dz/mem.splits_rad_z[i]) <= 1.0;
	      }
	    else
	      {
		inside=
		  pow(dx/mem.splits_rad_x[i],2) +
		  pow(dy/mem.splits_rad_y[i],2) +
		  pow(dz/mem.splits_rad_z[i],2) <= 1.0;
	      }
	  }
	if(inside)
	  many_split=max(many_split,mem.splits_many[i]);
      }
    if(mem.masks_init == 0)
      level_mask=mem.level_max;
    for(int i=0;i<mem.masks_init;i++)
      {
	bool inside=false;
	if(mem.masks_level_init[i] > level_mask)
	  {
	    double dx=x0-mem.masks_center_x_init[i];
	    double dy=y0-mem.masks_center_y_init[i];
	    double dz=z0-mem.masks_center_z_init[i];
	    if(mem.periodic)
	      {
		dx=min(min(abs(dx),abs(dx+1.0)),abs(dx-1.0));
		dy=min(min(abs(dy),abs(dy+1.0)),abs(dy-1.0));
		dz=min(min(abs(dz),abs(dz+1.0)),abs(dz-1.0));
	      }
	    if(mem.masks_square_init[i])
	      {
		inside=
		  abs(dx/mem.masks_rad_x_init[i]) <= 1.0 &&
		  abs(dy/mem.masks_rad_y_init[i]) <= 1.0 &&
		  abs(dz/mem.masks_rad_z_init[i]) <= 1.0;
	      }
	    else
	      {
		inside=
		  pow(dx/mem.masks_rad_x_init[i],2) +
		  pow(dy/mem.masks_rad_y_init[i],2) +
		  pow(dz/mem.masks_rad_z_init[i],2) <= 1.0;
	      }
	  }
	if(inside)
	  level_mask=max(level_mask,mem.masks_level_init[i]);
      }
    int many=min(many_split,Misc::pow(2,level_mask));
    if(split_to < 1) return many;
    if(!gen_part) 
      {
	count+=many*many*many;
	return many;
      }
    //
    if(many ==1)
      {
	try
	  {
	    par=new Particle(iphase,jfield);
	  }
	catch (bad_alloc& ba)
	  {
	    cerr << " single Particle allocation failure split particle ";
	    cerr << mem.p_mess->FractalRank << " " << count << " " << ba.what() << endl;
	    exit(0);
	  }
	frac.particle_list[count]=par;
	frac.particle_list[count]->set_pos(x0,y0,z0);
	frac.particle_list[count]->set_mass(m);
	count++;
      }
    else
      {
	double delta=1.0/(double)(many*mem.grid_length);
	double offset=-0.5*delta*(double)(many-1);
	double mm=m/pow((double)many,3);
	for(int iz=0;iz < many;iz++)
	  {
	    for(int iy=0;iy < many;iy++)
	      {
		for(int ix=0;ix < many;ix++)
		  {
		    try
		      {
			par=new Particle(iphase,jfield);
		      }
		    catch (bad_alloc& ba)
		      {
			cerr << " Multi Particle allocation failure split particle ";
			cerr << mem.p_mess->FractalRank << " " << count << " " << ba.what() << endl;
			exit(0);
		      }
		    frac.particle_list[count]=par;
		    frac.particle_list[count]->
		      set_pos(
			      x0+offset+delta*(double)ix,
			      y0+offset+delta*(double)iy,
			      z0+offset+delta*(double)iz);
		    frac.particle_list[count]->set_mass(mm);
		    count++;
		  }
	      }
	  }
      }
    return many;
  }
}
namespace FractalSpace
{
  template int split_particle(Fractal_Memory& mem,Fractal& frac,const double& x0,const double& y0,const double& z0,int& count,const double& m,const int& split_to,
const bool& gen_part);
}
