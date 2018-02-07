#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void scatter_particles(Fractal_Memory& mem,Fractal& frac)
  {
    assert(mem.FractalNodes==mem.p_mess->FractalNodes);
    frac.timing(-1,27);
    vector <double>timing(3);
    timing[0]=-mem.p_mess->Clock();
    ofstream& FF=frac.p_file->DUMPS;
    if(frac.get_periodic())
      {
	frac.wrap();
	add_pseudo_particles(mem,frac);
      }
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    int FractalNodes0=mem.FractalNodes0;
    int FractalNodes1=mem.FractalNodes1;
    int FractalNodes2=mem.FractalNodes2;
    vector <double> pos(3);
    vector <double> RealBox(6);
    vector <double> RealPBox(6);
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
    timing[0]+=mem.p_mess->Clock();
    timing[1]=-mem.p_mess->Clock();
    vector <int>posI(3);
    double SCALE=mem.grid_length;
    double DB=0.0;
    int DBI=0;
    if(mem.periodic)
      {
	DB=2.0/static_cast<double>(mem.grid_length);
	DBI=2;
      }
    // frac.particle_list.clear();
    clean_deque(frac.particle_list);
    int field_length=4;
    if(mem.calc_density_particle)
      field_length=5;
    if(mem.calc_shear)
      field_length=6;
    const bool change_it= field_length > 4;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=4;
    clean_deque(mem.p_mess->parts_tmp);
    int LOOPS=3;
    for(int LOOP=0;LOOP<LOOPS;LOOP++)
      {
	vector <int> counts_out(FractalNodes,0);
	vector <vector <int> > dataI_out(FractalNodes);
	vector <vector <double> > dataR_out(FractalNodes);
	vector <int> counts_in;
	vector <double> dataR_in;
	vector <int> dataI_in;
	for(int particle=0; particle < frac.get_number_particles_world(); ++particle)
	  {
	    if(particle % LOOPS == LOOP)
	      {
		Particle* P=frac.particle_list_world[particle];
		P->get_pos(pos);
		if(!(mem.periodic || vector_in_box(pos,mem.BigBox)))
		  continue;
		posI[0]=static_cast<int>((pos[0]+DB)*SCALE)-DBI;
		posI[1]=static_cast<int>((pos[1]+DB)*SCALE)-DBI;
		posI[2]=static_cast<int>((pos[2]+DB)*SCALE)-DBI;
		int FRZ=std::upper_bound(lowerZ.begin(),lowerZ.end(),posI[2])
		  -lowerZ.begin()-1;
		int FRY=std::upper_bound(lowerY[FRZ].begin(),lowerY[FRZ].end(),posI[1])
		  -lowerY[FRZ].begin()-1;
		int FRX=std::upper_bound(lowerX[FRZ][FRY].begin(),lowerX[FRZ][FRY].end(),posI[0])
		  -lowerX[FRZ][FRY].begin()-1;
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
	  }
	mem.p_mess->Send_Data_Some_How(5,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	dataR_out.clear();
	dataI_out.clear();
	int particle=0;
	Particle* pt=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		if(c == 0)
		  {
		    try
		      {
			pt=new Particle[counts_in[FR]];
			mem.p_mess->parts_tmp.push_back(pt);
		      }
		    catch(bad_alloc& ba)
		      {
			cerr << " bad particle scatter a " << " " << ba.what() << " " << FR << " " << counts_in[FR] << " " << frac.particle_list.size() << endl;
			exit(0);
		      }
		  }
		Particle* P=&pt[c];
		frac.particle_list.push_back(P);
		int p4=particle*4;
		P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],dataI_in[particle],FR);
		if(change_it)
		  P->field_resize(field_length);
		particle++;
	      }
	  }
      }
    timing[1]+=mem.p_mess->Clock();
    timing[2]=-mem.p_mess->Clock();
    frac.set_number_particles(frac.particle_list.size());
    // mem.TouchWhichBoxes.clear();
    clean_vector(mem.TouchWhichBoxes);
    integers=0;
    how_manyI=-1;
    how_manyR=-1;
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
    // mem.p_mess->parts_tmpp.clear();
    clean_deque(mem.p_mess->parts_tmpp);
    //    LOOPS=1;
    LOOPS=9;
    for(int LOOP=0;LOOP<LOOPS;LOOP++)
      {
	vector <int> counts_out(FractalNodes,0);
	vector <vector <int> > dataI_out(FractalNodes);
	vector <vector <double> > dataR_out(FractalNodes);
	vector <int> counts_in;
	vector <double> dataR_in;
	vector <int> dataI_in;
	for(int particle=0; particle < frac.get_number_particles(); ++particle)
	  {
	    if(particle % LOOPS == LOOP)
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
	  }
	mem.p_mess->Send_Data_Some_How(6,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	dataR_out.clear();
	dataI_out.clear();
	Particle* ppt=0;
	int particle=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		if(c == 0)
		  {
		    try
		      {
			ppt=new Particle[counts_in[FR]];
			mem.p_mess->parts_tmpp.push_back(ppt);
		      }
		    catch(bad_alloc& ba)
		      {
			cerr << " bad particle scatter b " << " " << ba.what() << " " << FR << " " << counts_in[FR] << " " << frac.particle_list.size() << endl;
			exit(0);
		      }
		  }
		Particle* P=&ppt[c];
		frac.particle_list.push_back(P);
		int p4=particle*4;
		P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],-1,FR);
		if(change_it)
		  P->field_resize(field_length);
		particle++;
	      }
	  }
      }
    timing[2]+=mem.p_mess->Clock();
    fprintf(mem.p_file->PFTime," scatter particles %10.3E %10.3E %10.3E \n",
	    timing[0],timing[1],timing[2]);
    frac.set_number_particles(frac.particle_list.size());
    frac.timing(1,27);
  }
}
