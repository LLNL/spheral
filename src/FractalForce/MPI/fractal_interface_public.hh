#ifndef _InterFacePublic_Defined_
#define _InterFacePublic_Defined_
namespace FractalSpace
{
  void doFractalForce(Fractal_Memory* PFM);
  Fractal_Memory* fractal_memory_create();
  void fractal_memory_delete(Fractal_Memory* PFM);
  void fractal_memory_content_delete(Fractal_Memory* PFM);
  void fractal_create(Fractal_Memory* PFM);
  void fractal_delete(Fractal_Memory* PFM);
  void addParticles(Fractal_Memory* PFM,int first,int total,
		    vector <double>& xmin,vector <double>& xmax,
		    vector <double>& xpos,vector <double>& ypos,
		    vector <double>& zpos,vector <double>& masses);
  void getField(Fractal_Memory& PFM,int first,int last,double G,
		vector <double>& xmin,vector <double>& xmax,
		vector <double>& pot,vector <double>& fx,
		vector <double>& fy,vector <double>& fz);
  /*
  void Fractal_Memory::fractal_memory_setup()
  void Fractal_Memory::setNumberParticles(int NP)
  void Fractal_Memory::setFractalNodes(int FR0,int FR1,int FR2)
  void Fractal_Memory::setPeriodic(bool per)
  void Fractal_Memory::setDebug(bool Db)
  void Fractal_Memory::setGridLength(int GL)
  void Fractal_Memory::setPadding(int PA)
  void Fractal_Memory::setLevelMax(int LM)
  void Fractal_Memory::setMinimumNumber(int MN)
  void Fractal_Memory::setHypreIterations(int MHI)
  void Fractal_Memory::setHypreTolerance(double HT)
  void Fractal_Memory::setBaseDirectory(string BD)
  void Fractal_Memory::setRunIdentifier(string RI)
  */
}
#endif
