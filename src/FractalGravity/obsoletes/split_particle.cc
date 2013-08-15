#include "libs.hh"
#include "classes.hh"
#include "nbody.hh"
#include "headers.hh"
void split_particle(Nbody& nbody,const double& x0,const double& y0,const double& z0,
		    int& count,const double& m)
{
  int level_split=0;
  int level_mask=0;
  for(int i=0;i<nbody.splits;i++)
    {
      bool inside=false;
      if(nbody.splits_level[i] > level_split)
	{
	  double dx=x0-nbody.splits_center_x[i];
	  double dy=y0-nbody.splits_center_y[i];
	  double dz=z0-nbody.splits_center_z[i];
	  inside=false;
	  if(nbody.periodic)
	    {
	      dx=min(min(abs(dx),abs(dx+1.0)),abs(dx-1.0));
	      dy=min(min(abs(dy),abs(dy+1.0)),abs(dy-1.0));
	      dz=min(min(abs(dz),abs(dz+1.0)),abs(dz-1.0));
	    }
	  if(nbody.splits_square[i])
	    {
	      inside=
		abs(dx/nbody.splits_rad_x[i]) <= 1.0 &&
		abs(dy/nbody.splits_rad_y[i]) <= 1.0 &&
		abs(dz/nbody.splits_rad_z[i]) <= 1.0;
	    }
	  else
	    {
	      inside=
		pow(dx/nbody.splits_rad_x[i],2) +
		pow(dy/nbody.splits_rad_y[i],2) +
		pow(dz/nbody.splits_rad_z[i],2) <= 1.0;
	    }
	}
      if(inside)
	level_split=max(level_split,nbody.splits_level[i]);
    }
  if(nbody.masks == 0)
    level_mask=nbody.level_max;
  for(int i=0;i<nbody.masks;i++)
    {
      bool inside=false;
      if(nbody.masks_level[i] > level_mask)
	{
	  double dx=x0-nbody.masks_center_x[i];
	  double dy=y0-nbody.masks_center_y[i];
	  double dz=z0-nbody.masks_center_z[i];
	  if(nbody.periodic)
	    {
	      dx=min(min(abs(dx),abs(dx+1.0)),abs(dx-1.0));
	      dy=min(min(abs(dy),abs(dy+1.0)),abs(dy-1.0));
	      dz=min(min(abs(dz),abs(dz+1.0)),abs(dz-1.0));
	    }
	  if(nbody.masks_square[i])
	    {
	      inside=
		abs(dx/nbody.masks_rad_x[i]) <= 1.0 &&
		abs(dy/nbody.masks_rad_y[i]) <= 1.0 &&
		abs(dz/nbody.masks_rad_z[i]) <= 1.0;
		}
	  else
	    {
	      inside=
		pow(dx/nbody.masks_rad_x[i],2) +
		pow(dy/nbody.masks_rad_y[i],2) +
		pow(dz/nbody.masks_rad_z[i],2) <= 1.0;
		}
	}
      if(inside)
	level_mask=max(level_mask,nbody.masks_level[i]);
    }
  int lev=min(level_split,level_mask);
  //
  if(lev ==0)
    {
      nbody.pos_x[count]=x0;
      nbody.pos_y[count]=y0;
      nbody.pos_z[count]=z0;
      //      cout << "pos0a " << count << " " << nbody.pos_x[count]<< " " << nbody.pos_y[count]<< " " << nbody.pos_z[count] << " " << m << endl;
      //       nbody.vel_x[count]=0.0;
      //       nbody.vel_y[count]=0.0;
      //       nbody.vel_z[count]=0.0;
      nbody.particle_mass[count]=m;
      count++;
  }
  else
    {
      int many=Misc::pow(2,lev);
      double delta=1.0/(double)(many*nbody.grid_length);
      double offset=-0.5*delta*(double)(many-1);
      double mm=m/pow(8.0,lev);
      for(int iz=0;iz < many;iz++)
	{
	  for(int iy=0;iy < many;iy++)
	    {
	      for(int ix=0;ix < many;ix++)
		{
		  nbody.pos_x[count]=x0+offset+delta*(double)ix;
		  nbody.pos_y[count]=y0+offset+delta*(double)iy;
		  nbody.pos_z[count]=z0+offset+delta*(double)iz;
		  //		  cout << "pos0b " << count << " " << many << " " << ix << " " << iy << " " << iz << " " << endl;
		  //		    nbody.pos_x[count]<< " " << nbody.pos_y[count]<< " " << nbody.pos_z[count] << " " << mm << endl;
		  // 		  nbody.vel_x[count]=0.0;
		  // 		  nbody.vel_y[count]=0.0;
		  // 		  nbody.vel_z[count]=0.0;
		  nbody.particle_mass[count]=mm;
		  count++;
		}
	    }
	}
    }
}
