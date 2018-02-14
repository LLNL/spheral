#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void fractal_force_init(Fractal_Memory* PFM)
  {
    // Fractal_Memory& fractal_memory=*PFM;
    // File* p_file=PFM->p_file;
    // ofstream& FileFractal=p_file->DUMPS;
    // p_file->note(true," a fractal_memory ");
    // bool root=PFM->p_mess->FractalRank == 0;
    PFM->calc_FractalNodes();
    // p_file->note(true," b fractal_memory ");
    // vector < vector <int> >BoxA=PFM->Boxes;
    // int FractalNodes=PFM->FractalNodes;
    PFM->calc_Buffers_and_more();
    PFM->calc_RealBoxes();
    // p_file->note(true," d fractal_memory ");
    // BoxA=PFM->Boxes;
    // for(int FR=0;FR < FractalNodes;FR++)
    //   {
    // 	if(root)
    // 	  FileFractal << "Box fracDD " << BoxA[FR][0] << " " << BoxA[FR][1] << " " << BoxA[FR][2] << " " 
    // 		      << BoxA[FR][3] << " " << BoxA[FR][4] << " " << BoxA[FR][5] << "\n";
    //   }
  }
}
 
