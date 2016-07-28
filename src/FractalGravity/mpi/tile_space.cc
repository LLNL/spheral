#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void gen_parts(vector <Particle*>& list_part);
  void gen_helpBoxes(vector <Box*>& helpBoxes);
  void gen_overlaps(vector <Box*>& helpBoxes,list <Box*>& list_box_big);
  int Particle::number_particles=0;
}
int main()
{
  using namespace FractalSpace;
  vector <Particle*> list_part;
  gen_parts(list_part);
  //
  vector <Box*> helpBoxes;
  gen_helpBoxes(helpBoxes);
  //
  list <Box*> list_box_big;

  vector <double> cb;
  cb.assign(3,0.5);
  vector <double> wb;
  wb.assign(3,0.5);
  Box* p_box_big= new Box(cb,wb);
  list_box_big.push_back(p_box_big);
  //
  gen_overlaps(helpBoxes,list_box_big);
  //
  return 1;
}
namespace FractalSpace
{
  void gen_parts(vector <Particle*>& list_part)
  {
    double x0=0.5;
    double y0=0.35;
    double z0=0.55;
    double r0=0.3;
    double rand_max=(double)RAND_MAX;
    int npoints=100000;
    int random=456;
    srand(random);
    double twopi=8.0*atan(1.0);
    vector <double> pos(3);
    for (int point=0;point<npoints;point++)
      {
	double r=r0*Fractal::my_rand(rand_max);
	double phi=twopi*Fractal::my_rand(rand_max);
	double ctheta=2.0*Fractal::my_rand(rand_max)-1.0;
	double stheta=sqrt(abs(1.0-ctheta*ctheta));
	pos[0]=x0+r*stheta*cos(phi);
	pos[1]=y0+r*stheta*sin(phi);
	pos[2]=z0+r*ctheta;
	Particle* p_part=new Particle(3,1);
	list_part.push_back(p_part);
	p_part->set_pos(pos);
      }
  }
  void gen_helpBoxes(vector <Box*>& helpBoxes)
  {
    int cubes=10;
    int cubes3=cubes*cubes*cubes;
    double acubes=cubes;
    double cubesinv=1.0/acubes;
    helpBoxes.resize(cubes3);
    vector <double> WC;
    WC.assign(3,cubesinv);
    vector <double> CC;
    CC.resize(3);
    int nC=0;
    for (int nz=0;nz < cubes;nz++)
      {
	CC[2]=(double)nz*cubesinv;
	for (int ny=0;ny < cubes;ny++)
	  {
	    CC[1]=(double)ny*cubesinv;
	    for (int nx=0;nx < cubes;nx++)
	      {
		CC[0]=(double)nx*cubesinv;
		helpBoxes[nC]->setW(WC);
		helpBoxes[nC]->setC(CC);
		nC++;
	      }
	  }
      }
  }
  void gen_overlaps(vector <Box*>& helpBoxes,list <Box*>& list_box_big)
  {
    int cubes3=(int)helpBoxes.size();
    for(int cube=0;cube < cubes3;cube++)
      {
	helpBoxes[cube]->overlaps.clear();
	for(list <Box*>::const_iterator box_b_itr=list_box_big.begin();box_b_itr!=list_box_big.end();box_b_itr++)
	  {
	    Box* p_box_b=*box_b_itr;
	    Box& box_b=*p_box_b;
	    if(helpBoxes[cube]->overlap(box_b))
	      helpBoxes[cube]->overlaps.push_back(p_box_b);
	}
    }
 }
}
