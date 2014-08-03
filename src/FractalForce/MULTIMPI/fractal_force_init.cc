#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void fractal_force_init(Fractal_Memory* p_fractal_memory)
  {
    Fractal_Memory& fractal_memory=*p_fractal_memory;
    File* p_file=p_fractal_memory->p_file;
    ofstream& FileFractal=p_file->DUMPS;
    p_file->note(true," a fractal_memory ");
    bool root=fractal_memory.p_mess->FractalRank == 0;
    fractal_memory.calc_FractalNodes();
    p_file->note(true," b fractal_memory ");
    vector < vector <int> >BoxA=p_fractal_memory->Boxes;
    int FractalNodes=p_fractal_memory->FractalNodes;
    /*
    for(int FR=0;FR < FractalNodes;FR++)
      {
	if(root)
	  FileFractal << "Box fracBB " << BoxA[FR][0] << " " << BoxA[FR][1] << " " << BoxA[FR][2] << " " 
		      << BoxA[FR][3] << " " << BoxA[FR][4] << " " << BoxA[FR][5] << "\n";
      }
    */
    fractal_memory.calc_Buffers_and_more();
    /*
    p_file->note(true," c fractal_memory ");
    BoxA=p_fractal_memory->Boxes;
    for(int FR=0;FR < FractalNodes;FR++)
      {
	if(root)
	  FileFractal << "Box fracCC " << BoxA[FR][0] << " " << BoxA[FR][1] << " " << BoxA[FR][2] << " " 
		      << BoxA[FR][3] << " " << BoxA[FR][4] << " " << BoxA[FR][5] << "\n";
      }
    */
    fractal_memory.calc_RealBoxes();
    p_file->note(true," d fractal_memory ");
    BoxA=p_fractal_memory->Boxes;
    for(int FR=0;FR < FractalNodes;FR++)
      {
	if(root)
	  FileFractal << "Box fracDD " << BoxA[FR][0] << " " << BoxA[FR][1] << " " << BoxA[FR][2] << " " 
		      << BoxA[FR][3] << " " << BoxA[FR][4] << " " << BoxA[FR][5] << "\n";
      }
    //    p_fractal_memory->p_fractal->redo(p_fractal_memory);

    /*
    if(p_fractal_memory->MPIrun)
      {
	p_file->note(true," e fractal_memory ");
	fractal_memory.p_mess->MPIStartup();
	p_file->note(true," f fractal_memory ");
	int length=fractal_memory.grid_length;
	bool periodic=fractal_memory.periodic;
	fractal_memory.p_mess->FFTWStartup(length,periodic);
      }
    else
      {
	fractal_memory.p_mess->start_x=0;
	fractal_memory.p_mess->length_x=fractal_memory.grid_length;
      }
    */
  }
}
 
