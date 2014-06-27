#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void scatter_particles(Fractal_Memory& mem,Fractal& frac)
  {
    assert(mem.FractalNodes==mem.p_mess->FractalNodes);
    ofstream& FF=frac.p_file->DUMPS;
    FF << " entered scatter particles " << "\n";
    if(frac.get_periodic())
      {
	frac.wrap();
	add_pseudo_particles(mem,frac);
      }
    if(!mem.MPIrun)
      return;
    double time0=mem.p_mess->Clock();
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    int FractalNodes0=mem.FractalNodes0;
    int FractalNodes1=mem.FractalNodes1;
    int FractalNodes2=mem.FractalNodes2;
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
    FF << "BIGBOX " << mem.BigBox[0] << " " << mem.BigBox[1] << " " << mem.BigBox[2] << " " << mem.BigBox[3] << " " << mem.BigBox[4] << " " << mem.BigBox[5] << "\n";
    int FNodesXY=FractalNodes0*FractalNodes1;
    vector <int> lowerZ(FractalNodes2);
    vector < vector <int> > lowerY(FractalNodes2);
    vector < vector < vector <int> > > lowerX(FractalNodes2);
    for(int FRZ=0;FRZ < FractalNodes2;FRZ++)
      {
	lowerZ[FRZ]=mem.LeftCorners[FRZ*FNodesXY][2];
	lowerY[FRZ].resize(FractalNodes1);
	lowerX[FRZ].resize(FractalNodes1);
	for(int FRY=0;FRY < FractalNodes1;FRY++)
	  {
	    lowerY[FRZ][FRY]=mem.LeftCorners[FRZ*FNodesXY+FRY*FractalNodes0][1];
	    lowerX[FRZ][FRY].resize(FractalNodes0);
	    for(int FRX=0;FRX < FractalNodes0;FRX++)
	      {
		lowerX[FRZ][FRY][FRX]=mem.LeftCorners[FRZ*FNodesXY+FRY*FractalNodes0+FRX][0];
	      }
	  }
      }
    double time1=mem.p_mess->Clock();
    vector <int>posI(3);
    double SCALE=mem.grid_length;
    double DB=0.0;
    int DBI=0;
    if(mem.periodic)
      {
	DB=2.0/static_cast<double>(mem.grid_length);
	DBI=2;
      }
    for(int particle=0; particle < frac.get_number_particles_world(); ++particle)
      {
	Particle* P=frac.particle_list_world[particle];
	P->get_pos(pos);
	if(!(mem.periodic || vector_in_box(pos,mem.BigBox)))
	  continue;
	posI[0]=static_cast<int>((pos[0]+DB)*SCALE)-DBI;
	posI[1]=static_cast<int>((pos[1]+DB)*SCALE)-DBI;
	posI[2]=static_cast<int>((pos[2]+DB)*SCALE)-DBI;
	vector <int>::iterator itza=lowerZ.begin();
	vector <int>::iterator itzb=lowerZ.end();
	int FRZ=std::upper_bound(itza,itzb,posI[2])-itza-1;
	vector <int>::iterator itya=lowerY[FRZ].begin();
	vector <int>::iterator ityb=lowerY[FRZ].end();
	int FRY=std::upper_bound(itya,ityb,posI[1])-itya-1;
	vector <int>::iterator itxa=lowerX[FRZ][FRY].begin();
	vector <int>::iterator itxb=lowerX[FRZ][FRY].end();
	int FRX=std::upper_bound(itxa,itxb,posI[0])-itxa-1;
	int FR=FRX+(FRY+FRZ*FractalNodes1)*FractalNodes0;
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
    double time2=mem.p_mess->Clock();
    vector <int> counts_in(FractalNodes);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=4;
    FF << "send stuff to other nodes a " << "\n";
    double time3=mem.p_mess->Clock();
    FF << "send stuff to other nodes b " << "\n";
    mem.p_mess->Send_Data_Some_How(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    double time4=mem.p_mess->Clock();
    FF << "send stuff to other nodes c " << "\n";
    dataR_out.clear();
    dataI_out.clear();
    frac.particle_list.resize(how_manyI);
    Particle* particles_tmp;
    try
      {
	particles_tmp=new Particle[how_manyI];
      }    
    catch(bad_alloc& ba)
      {
	cerr << " bad particle scatter " << " " << how_manyI << " " << ba.what() << endl;
	exit(0);
      }
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
    double time5=mem.p_mess->Clock();
    dataR_in.clear();
    dataI_in.clear();
    //    really_clear(dataR_in);
    //    really_clear(dataI_in);
    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    counts_out.assign(FractalNodes,0);
    counts_in.assign(FractalNodes,0);
    mem.TouchWhichBoxes.clear();
    integers=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	if(FractalRank == FR)
	  continue;
	if(overlap_boxes(mem.Boxes[FractalRank],mem.PBoxes[FR]))
	  {
	    mem.TouchWhichBoxes.push_back(FR);
	  }
      }
    int TBsize=mem.TouchWhichBoxes.size();
    FF << " IBOX " << mem.RealIBoxes[FractalRank][0] << " " << mem.RealIBoxes[FractalRank][1] << " " << mem.RealIBoxes[FractalRank][2] << " ";
    FF << mem.RealIBoxes[FractalRank][3] << " " << mem.RealIBoxes[FractalRank][4] << " " << mem.RealIBoxes[FractalRank][5] << "\n";
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	P->get_pos(pos);
	if(vector_in_box(pos,mem.RealIBoxes[FractalRank]))
	  continue;
	double pm=P->get_mass();
	for(int TB=0;TB<TBsize;TB++)
	  {
	    int FR=mem.TouchWhichBoxes[TB];
	    if(!vector_in_box(pos,mem.RealPBoxes[FR]))
	      continue;
	    dataR_out[FR].push_back(pos[0]);
	    dataR_out[FR].push_back(pos[1]);
	    dataR_out[FR].push_back(pos[2]);
	    dataR_out[FR].push_back(pm);
	    counts_out[FR]++;
	  }
      }
    double time6=mem.p_mess->Clock();
    FF << "send stuff to other nodes d " << "\n";
    double time7=mem.p_mess->Clock();
    FF << "send stuff to other nodes e " << "\n";
    mem.p_mess->Send_Data_Some_How(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    double time8=mem.p_mess->Clock();
    FF << "send stuff to other nodes f " << "\n";
    dataR_out.clear();
    dataI_out.clear();
    //    really_clear(dataR_out);
    //    really_clear(dataI_out);

    Particle* particles_tmpp=new Particle[(how_manyR/4)];
    mem.p_mess->parts_tmpp=particles_tmpp;
    particle=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    P=&particles_tmpp[particle];
	    assert(P);
	    frac.particle_list.push_back(P);
	    p4=particle*4;
	    P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],-1,FR);
	    if(change_it)
	      P->field_resize(field_length);
	    particle++;
	  }
      }
    frac.set_number_particles(frac.particle_list.size());
    double time9=mem.p_mess->Clock();
    fprintf(mem.p_file->PFTime," scatter particles %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E \n",
	    time1-time0,time2-time1,time3-time2,time4-time3,time5-time4,time6-time5,time7-time6,time8-time7,time9-time8,time9-time0);
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
