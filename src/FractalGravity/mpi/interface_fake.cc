// This is a fake interface to fractal_gravity for informational purposes
// so you can see what it needs. This is the simplest possible, you give
// fractal_gravity the positions and masses for all particle and it returns
// potential and forces and then forgets everything and releases all memory.
//
#include "libs.hh"
#include "classes.hh"
#include "jmofake.hh"
#include "headers.hh"
// These are all my header files, may have to be fixed up for consistency
//
// All my stuff sits it namespace FractalSpace
//
// I am faking your stuff just to show what I need
//
namespace FractalSpace
{
  using Spheral::jmodata;
  template <class N> void interface_fake(N& jmodata)
  {
    vector <double> pos_jmo(3,0.0);
    vector <double> pos_jvv(3,0.0);
    double scale_jmo=1.0;
    double G_jmo=1.0;
    jmodata.get_conversion(pos_jmo,scale_jmo,G_jmo);
    //This constructor copies fractal_memory info into fractal.
    //
    // The points (nodes) must fit into the unit box, so you must furnish a vector and a scalar so that
    // x_jvv=(x_jmo-pos_jmo0[0])/scale_jmo
    // y_jvv=(y_jmo-pos_jmo0[1])/scale_jmo
    // z_jvv=(z_jmo-pos_jmo0[2])/scale_jmo
    //
    // If you have isolated boundary conditions you can use your own masses
    // so "mass_scale_factor=1.0"
    // If you have periodic boundary conditions, the total mass must equal 3*Omega_0/(8*pi) 
    // Omega_0 is the value at the (expansion factor=1) epoch.
    // so "mass_scale_factor=total_mass_jmo/(3*Omega_0/(8*pi))"
    //
    // You need to provide a gravitational constant G_jmo as G_jvv=1.0
    //
    Fractal_Memory* p_fractal_memory= new (nothrow) Fractal_Memory;
    assert(p_fractal_memory);
    Fractal_Memory& fractal_memory=*p_fractal_memory; 
    //
    fractal_memory_parameters(fractal_memory);
    fractal_memory.number_particles=jmodata.number_particles;
    //
    Fractal* p_fractal=new (nothrow) Fractal(fractal_memory);
    assert(p_fractal);
    Fractal& fractal=*p_fractal;
    fractal.particle_list.resize(fractal_memory.number_particles);
    fix_memory(fractal,3,4);
    // This generates a vector of Particle pointers, each Particle has a phasespace vector of length 3(x,y,z) and a
    // field vector of length 4(pot,fx,fy,fz)=(0,0,0,0)
    double mass_jmo=0.0;
    double mass_jvv=0.0;
    int p=0;
    for (list <Node*>::const_itr node_itr=jmodata.list_nodes.begin();node_itr != jmodata.list_nodes.end();node_itr++)
      {
	Node& node=**node_itr;
	node.get_pos(pos_jmo);
	pos_jvv[0]=(pos_jmo[0]-pos_jmo0[0])/scale_jmo;
	pos_jvv[1]=(pos_jmo[1]-pos_jmo0[1])/scale_jmo;
	pos_jvv[2]=(pos_jmo[2]-pos_jmo0[2])/scale_jmo;
	mass_jmo=node.get_mass();
	mass_jvv=mass_jmo/mass_scale_factor;
	Particle& particle=*fractal.particle_list[p];
	particle.set_pos(pos_jvv);
	particle.set_mass(mass_jvv);
	p++;
      }
    fractal.timing(-2,0);
    fractal.timing(-1,29);
    fractal_gravity(fractal,fractal_memory);
    fractal.timing(1,29);
    fractal.timing(0,0);
    //
    vector <double> field_jvv(4);
    vector <double> field_jmo(4);
    double scale_jmo2=scale_jmo*scale_jmo;
    int p=0;
    for (list <Node*>::const_itr node_itr=jmodata.list_nodes.begin();node_itr != jmodata.list_nodes.end();node_itr++)
      {
	Node& node=**node_itr;
	Particle& particle=*fractal.particle_list[p];
	particle.get_field_pf(field_jvv);
	field_jmo[0]=field_jvv[0]*G_jmo/scale_jmo;
	field_jmo[1]=field_jvv[1]*G_jmo/scale_jmo2;
	field_jmo[2]=field_jvv[2]*G_jmo/scale_jmo2;
	field_jmo[3]=field_jvv[3]*G_jmo/scale_jmo2;
	node.set_field(field_jmo); 
	p++;
      }
    delete p_fractal;
    p_fractal=0;
    delete p_fractal_memory;
    p_fractal_memory=0;
  }
}
