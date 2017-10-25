namespace FractalSpace
{
  void fractal_create(Fractal_Memory* PFM);
  void fractal_delete(Fractal_Memory* PFM);
  template <class YDa,class YDb> void fractal_depopulate(Fractal_Memory* PFM, 
							 YDa& YourDataa, 
							 YDb& YourDatab);
  Fractal_Memory* fractal_interface_setup(bool periodic,
					  int grid_length,
					  int level_max,
					  int minimum_number,
					  int padding,
					  int tolHypre,
					  int maxitsHypre,
					  bool fractalDebug,
					  int FractalNodes0,
					  int FractalNodes1,
					  int FractalNodes2,
					  string BaseDirectory,
					  string RunIdentifier);
  void fractal_memory_content_delete(Fractal_Memory* PFM);
  void fractal_memory_interface_parameters(Fractal_Memory& mem,
					   bool periodic,
					   int grid_length,
					   int level_max,
					   int minimum_number,
					   int padding,
					   int tolHypre,
					   int maxitsHypre,
					   bool fractalDebug,
					   int FractalNodes0,
					   int FractalNodes1,
					   int FractalNodes2,
					   string BaseDirectory,
					   string RunIdentifier);
  template <class YDa,class YDb> void fractal_populate(Fractal_Memory* PFM, 
						       const int number_particles,
						       YDa& YourDataa, 
						       YDb& YourDatab);
}
