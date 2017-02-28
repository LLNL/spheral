#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void node_split(Fractal_Memory& mem,Fractal& frac)
  {
    if(mem.periodic)
      frac.wrap();
    if(mem.FractalNodes==1)
      return;
    int FractalNodes=mem.FractalNodes;
    int FractalNodes0=mem.FractalNodes0;
    int FractalNodes1=mem.FractalNodes1;
    int FractalNodes2=mem.FractalNodes2;
    int number_split=mem.number_split;
    int length=mem.grid_length;
    double alength=length;
    vector <double>xmin(3);
    xmin.assign(3,0.0);
    vector <double>xmax(3);
    xmax.assign(3,1.0);
    vector <double>pos(3);
    vector <double>tmp_x0;
    vector <double>tmp_y0;
    vector <double>tmp_z0;
    int nodeparts=frac.get_number_particles();
    int counts=0;
    for(int ni=0;ni<nodeparts;ni++)
      {
	frac.particle_list[ni]->get_pos(pos);
	if(mem.periodic || (
	  pos[0] >= xmin[0] && pos[0] <= xmax[0] &&
	  pos[1] >= xmin[1] && pos[1] <= xmax[1] &&
	  pos[2] >= xmin[2] && pos[2] <= xmax[2]))
	  {
	    tmp_x0.push_back(pos[0]);
	    tmp_y0.push_back(pos[1]);
	    tmp_z0.push_back(pos[2]);
	    counts++;
	  }
      }
    int nodeparts0=counts;
    double rand_max=(double)RAND_MAX;
    int particles[1]={nodeparts0};
    mem.p_mess->Find_Sum_INT(particles,1);
    int total_particles=particles[0];
    double aFN=FractalNodes;
    int sends=(nodeparts0*number_split)/(FractalNodes*total_particles);
    double anodeparts0=nodeparts0;
    int* counts_out=new int[FractalNodes];
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int ni=0;ni<sends;ni++)
	  {
	    int p=anodeparts0*Fractal::my_rand(rand_max);
	    dataR_out[FR].push_back(tmp_x0[p]);
	  }
	counts_out[FR]=sends;
      }
    int* counts_in=new int[FractalNodes];
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=0;
    int doubles=1;
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    mem.p_mess->Send_Data_Somewhere(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    dataR_out.clear();
    dataI_out.clear();
    sort(dataR_in.begin(),dataR_in.end());
    vector <double>xdiv(FractalNodes0+1);
    vector < vector <double> > ydiv;
    ydiv.resize(FractalNodes0+1);
    vector < vector < vector <double> > > zdiv;
    zdiv.resize(FractalNodes0+1);
    for(int FR0=0;FR0<FractalNodes0;FR0++)
      xdiv[FR0]=dataR_in[(FR0*how_manyR)/FractalNodes0]/aFN;
    xdiv[FractalNodes0]=dataR_in[how_manyR-1]/aFN;
    mem.p_mess->Find_Sum_DOUBLE(xdiv,FractalNodes0+1);
    delete [] counts_out;
    delete [] counts_in;
    //
    vector <double>tmp_x1;
    vector <double>tmp_y1;
    vector <double>tmp_z1;
    for(int FR0=0;FR0<FractalNodes0;FR0++)
      {
	int nodeparts1=0;
	for(int p=0;p<nodeparts0;p++)
	  {
	    if(tmp_x0[p] >= xdiv[FR0] && tmp_x0[p] < xdiv[FR0+1])
	      {
		tmp_x1.push_back(tmp_x0[p]);
		tmp_y1.push_back(tmp_y0[p]);
		tmp_z1.push_back(tmp_z0[p]);
		nodeparts1++;
	      }
	  }
	particles[1]=nodeparts1;
	mem.p_mess->Find_Sum_INT(particles,1);
	total_particles=particles[0];
	int sends=(nodeparts1*number_split)/(FractalNodes*total_particles);
	double anodeparts1=nodeparts1;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int ni=0;ni<sends;ni++)
	      {
		int p=anodeparts1*Fractal::my_rand(rand_max);
		dataR_out[FR].push_back(tmp_y1[p]);
	      }
	    counts_out[FR]=sends;
	  }
	int* counts_in=new int[FractalNodes];
	vector <double> dataR_in;
	vector <int> dataI_in;
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=0;
	int doubles=1;
	mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
	mem.p_mess->Send_Data_Somewhere(counts_out,counts_in,integers,doubles,
					dataI_out,dataI_in,how_manyI,
					dataR_out,dataR_in,how_manyR);
	dataR_out.clear();
	dataI_out.clear();
	sortit.assign(dataR_in,dataR_in+how_manyR);
	sort(sortit.begin(),sortit.end());
	ydiv[FR0].resize(FractalNodes1+1);
	zdiv[FR0].resize(FractalNodes1+1);
	for(int FR1=0;FR1<FractalNodes1;FR1++)
	  ydiv[FR0][FR1]=sortit[(FR1*how_manyR)/FractalNodes1]/aFN;
	ydiv[FR0][FractalNodes1]=sortit[how_manyR-1]/aFN;
	mem.p_mess->Find_Sum_DOUBLE(ydiv[FR0],FractalNodes0+1);
	delete [] counts_out;
	delete [] counts_in;
	delete [] dataI_in;
	delete [] dataR_in;
	vector <double>tmp_x2;
	vector <double>tmp_y2;
	vector <double>tmp_z2;
	for(int FR1=0;FR1<FractalNodes1;FR1++)
	  {
	    int nodeparts2=0;
	    for(int p=0;p<nodeparts1;p++)
	      {
		if(tmp_y1[p] >= ydiv[FR0][FR1] && tmp_y1[p] < ydiv[FR0][FR1+1])
		  {
		    tmp_x2[nodeparts2]=tmp_x1[p];
		    tmp_y2[nodeparts2]=tmp_y1[p];
		    tmp_z2[nodeparts2]=tmp_z1[p];
		    nodeparts2++;
		  }
	      }
	    particles[1]=nodeparts2;
	    mem.p_mess->Find_Sum_INT(particles,1);
	    total_particles=particles[0];
	    int sends=(nodeparts2*number_split)/(FractalNodes*total_particles);
	    double anodeparts2=nodeparts2;
	    zdiv[FR0][FR1].resize(FractalNodes);
	    for(int FR=0;FR<FractalNodes;FR++)
	      {
		for(int ni=0;ni<sends;ni++)
		  {
		    int p=anodeparts2*Fractal::my_rand(rand_max);
		    dataR_out[FR].push_back(tmp_z1[p]);
		  }
		counts_out[FR]=sends;
	      }
	    int* counts_in=new int[FractalNodes];
	    vector <double> dataR_in;
	    vector <int> dataI_in;
	    int how_manyI=-1;
	    int how_manyR=-1;
	    int integers=0;
	    int doubles=1;
	    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
	    mem.p_mess->Send_Data_Somewhere(counts_out,counts_in,integers,doubles,
					    dataI_out,dataI_in,how_manyI,
					    dataR_out,dataR_in,how_manyR);
	    dataR_out.clear();
	    dataI_out.clear();
	    sortit.assign(dataR_in,dataR_in+how_manyR);
	    sort(sortit.begin(),sortit.end());
	    zdiv[FR0][FR1].resize(FractalNodes2+1);
	    for(int FR2=0;FR2<FractalNodes2;FR2++)
	      zdiv[FR0][FR1][FR2]=sortit[(FR2*how_manyR)/FractalNodes2]/aFN;
	    zdiv[FR0][FR1][FractalNodes2]=sortit[how_manyR-1]/aFN;
	    mem.p_mess->Find_Sum_DOUBLE(zdiv[FR0][FR1],FractalNodes1+1);
	    delete [] counts_out;
	    delete [] counts_in;
	    delete [] dataI_in;
	    delete [] dataR_in;
	  }
      }
    for(int FR0=0;FR0<FractalNodes0;FR0++)
      {
	int n0a=xdiv[FR0]*alength;
	if(FR0 == 0)
	  n0a=0;
	int n0b=static_cast<int>(xdiv[FR0+1]*alength)-1;
	if(FR0==FractalNodes0-1)
	  n0b=length-1;
	for(int FR1=0;FR1<FractalNodes1;FR1++)
	  {
	    int n1a=ydiv[FR0][FR1]*alength;
	    if(FR1 ==0)
	      n1a=0;
	    int n1b=static_cast<int>(ydiv[FR0][FR1+1]*alength)-1;
	    if(FR1 == FractalNodes1-1)
	      n1b=length-1;
	    for(int FR2=0;FR2<FractalNodes2;FR2++)
	      {
		int n2a=zdiv[FR0][FR1][FR2]*alength;
		if(FR2 == 0)
		  n2a=0;
		int n2b=static_cast<int>(zdiv[FR0][FR1][FR2+1]*alength)-1;
		if(FR2 == FractalNodes2-1)
		  n2b=length-1;
		int FR=FR0+(FR1+FR2*FractalNodes1)*FractalNodes0;
		mem.Boxes[FR][0]=n0a;
		mem.Boxes[FR][1]=n0b;
		mem.Boxes[FR][2]=n1a;
		mem.Boxes[FR][3]=n1b;
		mem.Boxes[FR][4]=n2a;
		mem.Boxes[FR][5]=n2b;
	      }
	  }
      }
    mem.calc_Buffers_and_more();
    mem.calc_RealBoxes();
    frac.redo(mem);
  }
}
