#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void scatter_particles(Fractal_Memory& mem,Fractal& frac)
  {
    assert(mem.FractalNodes==mem.p_mess->FractalNodes);
    ofstream& FF=frac.p_file->FileFractal;
    frac.p_file->FileFractal << " entered scatter particles " << endl;
    if(frac.get_periodic())
      {
	frac.wrap();
	add_pseudo_particles(mem,frac);
      }
    if(!mem.MPIrun)
      return;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <double> pos_left(3);
    vector <double> pos_right(3);
    vector <double> pos_left_safe(3);
    vector <double> pos_right_safe(3);
    vector <double> pos(3);
    vector <double> RealBox(6);
    vector <double> RealPBox(6);
    vector <int> counts_out(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    left_right(frac,pos_left,pos_right);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      frac.particle_list[particle]->set_world(FractalRank,particle);
    frac.particle_list_world=frac.particle_list;
    frac.set_number_particles_world(frac.get_number_particles());
    FF << " safe " ;
    for(int ni=0;ni<3;ni++)
      {
	double delta=(pos_right[ni]-pos_left[ni])*0.5*(1.0+0.001);
	double mid=(pos_right[ni]+pos_left[ni])*0.5;
	pos_right_safe[ni]=mid+delta;
	pos_left_safe[ni]=mid-delta;
	FF << " " << pos_left[ni] << " " << pos_right[ni];
	FF << " " << pos_left_safe[ni] << " " << pos_right_safe[ni];
      }
    FF << endl;
    int blocks=23;
    double ablocks=blocks;
    ablocks=1.0/ablocks;
    int blocks3=blocks*blocks*blocks;
    vector < vector <double> > BigBlocks(blocks3);
    int ni=0;
    for(int ni2=0;ni2<blocks;ni2++)
      {
	double a2=static_cast<double>(ni2)*ablocks;
	for(int ni1=0;ni1<blocks;ni1++)
	  {
	    double a1=static_cast<double>(ni1)*ablocks;
	    for(int ni0=0;ni0<blocks;ni0++)
	      {
		double a0=static_cast<double>(ni0)*ablocks;
		BigBlocks[ni].resize(6);
		BigBlocks[ni][0]=a0;
		BigBlocks[ni][2]=a1;
		BigBlocks[ni][4]=a2;
		BigBlocks[ni][1]=BigBlocks[ni][0]+ablocks;;
		BigBlocks[ni][3]=BigBlocks[ni][2]+ablocks;;
		BigBlocks[ni][5]=BigBlocks[ni][4]+ablocks;;
		ni++;
	      }
	  }
      }
    double dx=pos_right_safe[0]-pos_left_safe[0];
    double dy=pos_right_safe[1]-pos_left_safe[1];
    double dz=pos_right_safe[2]-pos_left_safe[2];
    for(int ni=0;ni<blocks3;ni++)
      {
	BigBlocks[ni][0]=pos_left_safe[0]+BigBlocks[ni][0]*dx;
	BigBlocks[ni][1]=pos_left_safe[0]+BigBlocks[ni][1]*dx;
	BigBlocks[ni][2]=pos_left_safe[1]+BigBlocks[ni][2]*dy;
	BigBlocks[ni][3]=pos_left_safe[1]+BigBlocks[ni][3]*dy;
	BigBlocks[ni][4]=pos_left_safe[2]+BigBlocks[ni][4]*dz;
	BigBlocks[ni][5]=pos_left_safe[2]+BigBlocks[ni][5]*dz;
      }
    vector < vector <int> > LookHere(blocks3);
    for(int bl=0;bl<blocks3;bl++)
      {
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    if(overlap_boxes(mem.RealPBoxes[FR],BigBlocks[bl]))
	      {
		LookHere[bl].push_back(FR);
	      }
	  }
      }
    FF << " finished with big blocks" << endl;
    double bblocks=blocks;
    double dxinv=bblocks/dx;		
    double dyinv=bblocks/dy;		
    double dzinv=bblocks/dz;
    for(int particle=0; particle < frac.get_number_particles_world(); ++particle)
      {
	//	bool success=false;
	Particle* P=frac.particle_list_world[particle];
	P->get_pos(pos);
	int ni0=(pos[0]-pos_left_safe[0])*dxinv;
	int ni1=(pos[1]-pos_left_safe[1])*dyinv;
	int ni2=(pos[2]-pos_left_safe[2])*dzinv;
	int ni=ni0+(ni1+ni2*blocks)*blocks;
	int looks=LookHere[ni].size();
	for(int look=0;look<looks;look++)
	  {
	    int FR=LookHere[ni][look];
	    if(!vector_in_box(pos,mem.RealPBoxes[FR]))
	      continue;
	    int part=particle;
	    if(!vector_in_box(pos,mem.RealBoxes[FR]))
	      part=-particle-1;
	    dataR_out[FR].push_back(pos[0]);
	    dataR_out[FR].push_back(pos[1]);
	    dataR_out[FR].push_back(pos[2]);
	    dataR_out[FR].push_back(P->get_mass());
	    dataI_out[FR].push_back(part);
	    counts_out[FR]++;
	  }
      }

    vector <int> counts_in(FractalNodes);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=4;
    frac.timing(-1,33);
    mem.p_mess->Full_Stop();
    frac.timing(1,33);
    FF << "send stuff to other nodes a " << endl;
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    FF << "send stuff to other nodes b " << endl;
    mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    FF << "send stuff to other nodes c " << endl;
    dataR_out.clear();
    dataI_out.clear();
    frac.particle_list.resize(how_manyI);
    Particle* particles_tmp=new Particle[how_manyI];
    mem.p_mess->parts_tmp=particles_tmp;
    int field_length=4;
    if(mem.calc_density_particle)
      field_length=5;
    if(mem.calc_shear)
      field_length=6;
    const bool change_it= field_length > 4;
    vector <double>testRBox(6);
    //    testRBox=mem.RealPBoxes[FractalRank];
    //    vector <double>testRPos(3);
    int particle=0;
    int p4=-1;
    Particle* P=0;
    //    bool allok=true;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	FF << " testing " << FR << endl;
	for(int c=0;c<counts_in[FR];c++)
	  {
	    //	    FF << " testing " << FR << " " << c << endl;
	    P=&particles_tmp[particle];
	    assert(P);
	    frac.particle_list[particle]=P;
	    p4=particle*4;
	    P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],dataI_in[particle],FR);
	    //
	    /*
	    testRPos[0]=dataR_in[p4];
	    testRPos[1]=dataR_in[p4+1];
	    testRPos[2]=dataR_in[p4+2];
	    if(!vector_in_box(testRPos,testRBox))
	      {
		allok=false;
		FF << " Outside the Box " << FR << " " << c << " " << testRPos[0] << " " << testRPos[1] << " " << testRPos[2];
		FF << " " << testRBox[0] << " " << testRBox[1] << " " << testRBox[2];
		FF << " " << testRBox[3] << " " << testRBox[4] << " " << testRBox[5] << endl;
	      }
	    */
	    //
	    if(change_it)
	      P->field_resize(field_length);
	    particle++;
	  }
      }
    frac.set_number_particles(particle);
    //
    //    mem.p_mess->Full_Stop();
    //    assert(0);
  }
}
