#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  void doFractalForce(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    PF->set_steps(PFM->steps);
    PF->timing(-2,0);
    PF->timing(-1,49);
    fractal_force(*PF,*PFM);
    PFM->steps++;
    PF->timing(1,49);
    PF->timing(0,0);
    PF->timing_lev(0,0);
  }
  Fractal_Memory* fractal_memory_create()
  {
    Fractal_Memory* PFM= new Fractal_Memory;
    return PFM;
  }
  void fractal_memory_delete(Fractal_Memory* PFM)
  {
    delete PFM;
    PFM=0;
  }
  void Fractal_Memory::fractal_memory_setup()
  {
    MPIrun=FractalNodes > 1;
    global_level_max=level_max;
    /***********************/
    //Leave these parameters alone
    min_hypre_group_size=1;
    //    min_hypre_group_size=-1;
    // minum group size to use hypre
    new_points_gen=9;
    //Generate this many Points in each go
    steps=0;
    momentum_conserve=false;
    amnesia=true; // (true) forget everything after you are done. (false) remember everything.
    mind_wipe=false; // (true) delete everything and then come back without calculating anything.
    fixed_potential=false; // (true) use the fixed potential.
    calc_shear=false;// (true) if we calculate shear of force field
    calc_density_particle=false;
    do_vel=false;
    start_up=false;
    halo_fixed=false;
    //
    splits=0;
    masks=0;
    //
    masks_init=0;
    //
    p_file->note(true," a fractal_memory ");
    calc_FractalNodes();
    p_file->note(true," b fractal_memory ");
    calc_Buffers_and_more();
    p_file->note(true," c fractal_memory ");
    calc_RealBoxes();
    p_file->note(true," d fractal_memory ");
  }
  void fractal_memory_content_delete(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    Particle* P=PFM->p_mess->parts_interface;
    delete [] P;
    P=0;
    delete PF;
    PF=0;
    Mess* p_mess=PFM->p_mess;
    delete p_mess;
    p_mess=0;
    File* p_file=PFM->p_file;
    delete p_file;
    p_file=0;
  }
  void fractal_create(Fractal_Memory* PFM)
  {
    Fractal* PF=new Fractal(*PFM);
    PFM->p_fractal=PF;
    Particle* Parts_in=new Particle[PFM->number_particles];
    PFM->p_mess->parts_interface=Parts_in;
  }
  void addParticles(Fractal_Memory* PFM,int first,int total,
		    vector <double>& xmin,vector <double>& xmax,
		    vector <double>& xpos,vector <double>& ypos,
		    vector <double>& zpos,vector <double>& masses)
  {
    vector <double> pos(3);
    double dinv=1.0/(xmax[0]-xmin[0]);
    for(int ni=0;ni<total;ni++)
      {
	pos[0]=(xpos[ni]-xmin[0])*dinv;
	pos[1]=(ypos[ni]-xmin[1])*dinv;
	pos[2]=(zpos[ni]-xmin[2])*dinv;
	Particle* P=&PFM->p_mess->parts_interface[ni+first];
	PFM->p_fractal->particle_list.push_back(P);
	P->set_pos(pos);
	P->set_mass(masses[ni]);
      }
  }
  void getField(Fractal_Memory* PFM,int first,int total,double G,
		vector <double>& xmin,vector <double>& xmax,
		vector <double>& pot,vector <double>& fx,
		vector <double>& fy,vector <double>& fz)
  {
    
    vector <double> potforce(4);
    double dinv=1.0/(xmax[0]-xmin[0]);
    double convpot=G*dinv;
    double convforce=G*dinv*dinv;
    for(int ni=0;ni<total;ni++)
      {
	PFM->p_fractal->particle_list[ni+first]->get_field_pf(potforce);
	pot[ni]=potforce[0]*convpot;
	fx[ni]=potforce[1]*convforce;
	fy[ni]=potforce[2]*convforce;
	fz[ni]=potforce[3]*convforce;
      }
  }
  void fractal_delete(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    Particle* P=PFM->p_mess->parts_interface;
    delete [] P;
    P=0;
    delete PF;
    PF=0;
  }
  void Fractal_Memory::setNumberParticles(int NP)
  {
    number_particles=NP;
  }
  void Fractal_Memory::setFractalNodes(int FR0,int FR1,int FR2)
  {
    FractalNodes0=FR0;
    FractalNodes1=FR1;
    FractalNodes1=FR2;
    FractalNodes=FR0*FR1*FR2;
    MPIrun=FractalNodes > 1;
  }
  void Fractal_Memory::setPeriodic(bool per)
  {
    periodic=per;
  }
  void Fractal_Memory::setDebug(bool Db)
  {
    debug=Db;
  }
  void Fractal_Memory::setGridLength(int GL)
  {
    grid_length=GL;
  }
  void Fractal_Memory::setPadding(int PA)
  {
    PA=min(PA,1);
    PA=max(PA,-1);
    padding=PA;
  }
  void Fractal_Memory::setLevelMax(int LM)
  {
    level_max=LM;
  }
  void Fractal_Memory::setMinimumNumber(int MN)
  {
    minimum_number=MN;
  }
  void Fractal_Memory::setHypreIterations(int MHI)
  {
    maxits=MHI;
  }
  void Fractal_Memory::setHypreTolerance(double HT)
  {
    epsilon_sor=HT;
  }
  void Fractal_Memory::setBaseDirectory(string BD)
  {
    BaseDirectory=BD;
  }
  void Fractal_Memory::setRunIdentifier(string RI)
  {
    RUN=RI;
  }
}
