#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::vel="null";
  //
  void velocities(Fractal_Memory& mem,Fractal& frac)
  {
    ofstream& FileFractal= frac.p_file->DUMPS;
    //    ofstream& FileFractal= frac.p_file->FileFractal;
    FileFractal << "not yet " << &mem << " " << &frac << "\n";    
  }
}
