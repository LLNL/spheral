#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void assign_density(Group& group, Fractal& fractal)
  {
    ofstream& FileFractal=fractal.p_file->FileFractal;
    if(fractal.get_debug()) FileFractal << "enter assign density " << &group << " " << group.get_level() << endl;
    //--------------------------------------------------------------------------------------------------------------------------------
    // equivalength grid_length if total volume done at highest resolution
    //--------------------------------------------------------------------------------------------------------------------------------
    const  double scale=(double)(fractal.get_grid_length()*Misc::pow(2,fractal.get_level_max()));
    //--------------------------------------------------------------------------------------------------------------------------------
    //inverse width of a cell at current level, width at level=level_max is unity
    //--------------------------------------------------------------------------------------------------------------------------------
    const double d_inv=pow(2.0,group.get_level()-fractal.get_level_max());
    vector  <double> dens(8,0.0);
    double d_x,d_y,d_z;
    vector <double> pos(3);
    //--------------------------------------------------------------------------------------------------------------------------------
    // loop over all points in the group
    //--------------------------------------------------------------------------------------------------------------------------------
    //
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(point.list_particles.empty())
	  continue;
	//--------------------------------------------------------------------------------------------------------------------------------
	// loop over all particles associated with a point. only points with real_pointer=0,1,3,4,9,10,12,13
	// can have particles associated.
	//--------------------------------------------------------------------------------------------------------------------------------
	dens.assign(8,0.0);
	for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
	  {
	    Particle& particle=**particle_itr;
	    if(particle.get_mass() == 0.0)
	      continue;
	    //--------------------------------------------------------------------------------------------------------------------------------
	    //find particle postion in cell
	    //--------------------------------------------------------------------------------------------------------------------------------
	    particle.get_pos(pos);
	    point.get_deltas(pos,d_x,d_y,d_z,scale,d_inv);
	    //--------------------------------------------------------------------------------------------------------------------------------
	    // add mass to the eight corners
	    //--------------------------------------------------------------------------------------------------------------------------------
	    Misc::add_dens<double>(dens,particle.get_mass(),d_x,d_y,d_z);
	  }
	//--------------------------------------------------------------------------------------------------------------------------------
	// add mass to the eight points forming the corners. These eight points will always all exist
	// because of the restriction of which points can have particles
	//--------------------------------------------------------------------------------------------------------------------------------
	//	    FileFractal << " testing " << endl;
	//	    point.dump();
	point.add_density_at_points<double>(dens);
      }
    //--------------------------------------------------------------------------------------------------------------------------------
    // scale from mass at point to density at point
    //--------------------------------------------------------------------------------------------------------------------------------
    if(group.get_set_scaling())
      {
	//--------------------------------------------------------------------------------------------------------------------------------
	// scaling is the inverse volume of a cell
	//--------------------------------------------------------------------------------------------------------------------------------
	double scaling=(double)(Misc::pow(fractal.get_grid_length(),3))*pow(8.0,group.get_level());
	//      FileFractal << " scaling " << scaling << " " << group.get_level() << endl;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    //--------------------------------------------------------------------------------------------------------------------------------
	    // use scaling from mass per point to real density
	    //--------------------------------------------------------------------------------------------------------------------------------
	    (*point_itr)->scale_density_point(scaling);
	    //	    if((*point_itr)->get_density_point() != 0.0) 
	    //	      (*point_itr)->dumpd();
	  }
      }
    double d_0=fractal.get_density_0();
    //  FileFractal << "density zero " << d_0 << endl;
    //--------------------------------------------------------------------------------------------------------------------------------
    // subtract mean density
    //--------------------------------------------------------------------------------------------------------------------------------
    if(group.get_set_dens() && d_0 != 0.0)
      group.subtract_density(d_0);
    if(fractal.get_debug()) FileFractal << "leaving assign density" << endl;
    //
  }
}
