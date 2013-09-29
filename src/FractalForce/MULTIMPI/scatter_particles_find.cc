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
    vector <double> pos(3);
    vector <double> RealBox(6);
    vector <double> RealPBox(6);
    vector <int> counts_out(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      frac.particle_list[particle]->set_world(FractalRank,particle);
    frac.particle_list_world=frac.particle_list;
    frac.set_number_particles_world(frac.get_number_particles());
    FF << "BIGBOX " << mem.BigBox[0] << " " << mem.BigBox[1] << " " << mem.BigBox[2] << " " << mem.BigBox[3] << " " << mem.BigBox[4] << " " << mem.BigBox[5] << endl;
    for(int FR=0;FR<FractalNodes;FR++)
      FF << " LEFT " << FR << " " << mem.LeftCorners[FR][0]<< " " << mem.LeftCorners[FR][1] << " " << mem.LeftCorners[FR][2] << endl;
    vector <int>posI(3);
    double SCALE=mem.grid_length;
    double DB=0.0;
    double DBI=0;
    if(mem.periodic)
      {
	DB=1.0/static_cast<double>(mem.grid_length);
	DBI=1;
      }
    int FX=mem.FractalNodes0;
    int FXY=mem.FractalNodes0*mem.FractalNodes1;
    for(int particle=0; particle < frac.get_number_particles_world(); ++particle)
      {
	Particle* P=frac.particle_list_world[particle];
	P->get_pos(pos);
	//	FF << " pos " << particle << " " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
	if(!(mem.periodic || vector_in_box(pos,mem.BigBox)))
	  continue;
	posI[0]=static_cast<int>((pos[0]+DB)*SCALE)-DBI;
	posI[1]=static_cast<int>((pos[1]+DB)*SCALE)-DBI;
	posI[2]=static_cast<int>((pos[2]+DB)*SCALE)-DBI;
	//	FF << " posI " << particle << " " << posI[0] << " " << posI[1] << " " << posI[2] << endl;
	vector < vector <int> >::iterator ita=mem.LeftCorners.begin();
	vector < vector <int> >::iterator itb=mem.LeftCorners.end();
	vector < vector <int> >::iterator itFR=std::upper_bound(ita,itb,posI,compare_vectorsZ);
	advance(itFR,-FXY);
	ita=itFR;
	itb=itFR;
	advance(itb,FXY);
	itFR=std::upper_bound(ita,itb,posI,compare_vectorsY);
	advance(itFR,-FX);
	ita=itFR;
	itb=itFR;
	advance(itb,FX);
	itFR=std::upper_bound(ita,itb,posI,compare_vectorsX);
	advance(itFR,-1);
	int FR=itFR-mem.LeftCorners.begin();
	//	FF << " " << FR << " " << mem.RealBoxes[FR][0] << " " << mem.RealBoxes[FR][1] << " " << mem.RealBoxes[FR][2] << " " << mem.RealBoxes[FR][3] << " " << mem.RealBoxes[FR][4] << " " << mem.RealBoxes[FR][5] << endl;
	//	FF << " " << FR << " " << mem.Boxes[FR][0] << " " << mem.Boxes[FR][2] << " " << mem.Boxes[FR][4] << endl;
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
    vector <int> counts_in(FractalNodes);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=4;
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
    int particle=0;
    int p4=-1;
    Particle* P=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	FF << " testing a " << FR << endl;
	for(int c=0;c<counts_in[FR];c++)
	  {
	    P=&particles_tmp[particle];
	    assert(P);
	    frac.particle_list[particle]=P;
	    p4=particle*4;
	    P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],dataI_in[particle],FR);
	    if(change_it)
	      P->field_resize(field_length);
	    particle++;
	  }
      }
    frac.set_number_particles(particle);
    dataR_in.clear();
    dataI_in.clear();
    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    counts_out.assign(FractalNodes,0);
    counts_in.assign(FractalNodes,0);
    mem.TouchWhichBoxes.clear();
    for(int FR=0;FR<FractalNodes;FR++)
      {
	if(FractalRank == FR)
	  continue;
	if(overlap_boxes(mem.Boxes[FractalRank],mem.PBoxes[FR]))
	  {
	    mem.TouchWhichBoxes.push_back(FR);
	    //	    FF << " touch a " << FR << endl;
	  }
      }
    int TBsize=mem.TouchWhichBoxes.size();
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	P->get_pos(pos);
	//	FF << " neighs a" << particle << endl;
	if(vector_in_box(pos,mem.RealIBoxes[FractalRank]))
	  continue;
	//	FF << " neighs a" << particle << endl;
	int part=-particle-1;
	double pm=P->get_mass();
	for(int TB=0;TB<TBsize;TB++)
	  {
	    int FR=mem.TouchWhichBoxes[TB];
	    //	    FF << " neighs B " << TB << " " << FR << endl;
	    if(!vector_in_box(pos,mem.RealPBoxes[FR]))
	      continue;
	    dataR_out[FR].push_back(pos[0]);
	    dataR_out[FR].push_back(pos[1]);
	    dataR_out[FR].push_back(pos[2]);
	    dataR_out[FR].push_back(pm);
	    dataI_out[FR].push_back(part);
	    counts_out[FR]++;
	  }
      }
    FF << "send stuff to other nodes d " << endl;
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    FF << "send stuff to other nodes e " << endl;
    mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    FF << "send stuff to other nodes f " << endl;
    dataR_out.clear();
    dataI_out.clear();

    Particle* particles_tmpp=new Particle[how_manyI];
    mem.p_mess->parts_tmpp=particles_tmpp;
    particle=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	//	FF << " testing b " << FR << endl;
	for(int c=0;c<counts_in[FR];c++)
	  {
	    P=&particles_tmpp[particle];
	    assert(P);
	    frac.particle_list.push_back(P);
	    p4=particle*4;
	    P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],dataI_in[particle],FR);
	    if(change_it)
	      P->field_resize(field_length);
	    particle++;
	  }
      }
    frac.set_number_particles(frac.particle_list.size());
  }
  bool compare_vectorsZ(vector <int> veca,vector <int> vecb)
  {
    return veca[2]-vecb[2] < 0;
  }
  bool compare_vectorsY(vector <int> veca,vector <int> vecb)
  {
    return veca[1]-vecb[1] < 0;
  }
  bool compare_vectorsX(vector <int> veca,vector <int> vecb)
  {
    return veca[0]-vecb[0] < 0;
  }
}
