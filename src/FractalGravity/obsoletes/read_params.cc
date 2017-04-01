#include "libs.hh"
#include "classes.hh"
#include "nbody.hh"
#include "headers.hh"
void read_params(Nbody& nbody)
{
  ifstream InputFile;
  //  string FileNameIn;
  //  string FileNameOut;
  //  cin >> FileNameIn;
  InputFile.open("start.in");
  if(!InputFile)
    assert(0);
  int random_int;
  InputFile >> nbody. number_steps_total;
  InputFile >> nbody. number_steps_out;
  InputFile >> nbody. number_particles;
  InputFile >> nbody. grid_length;
  InputFile >> nbody. moat_0;
  InputFile >> nbody. periodic;
  InputFile >> nbody. random_initial;
  InputFile >> nbody. minimum_number;
  InputFile >> nbody. level_max;
  InputFile >> nbody. tweaks;
  InputFile >> nbody. padding;
  InputFile >> nbody. off_x;
  InputFile >> nbody. off_y;
  InputFile >> nbody. off_z;
  InputFile >> nbody. force_max;
  InputFile >> nbody. omega_0;
  InputFile >> nbody. omega_lambda;
  InputFile >> nbody. omega_b;
  InputFile >> nbody. h;
  InputFile >> nbody. box_length;
  InputFile >> nbody. norm_what;
  InputFile >> nbody. spectrum_number;
  InputFile >> nbody. highest_level_init;
  InputFile >> nbody. power_slope;
  InputFile >> nbody. norm_scale;
  InputFile >> nbody. sigma_init;
  InputFile >> nbody. step_length;
  InputFile >> nbody. pexp;
  InputFile >> nbody. redshift_start;
  InputFile >> random_int;
  srand(random_int);
  /*
  InputFile >> FileNameOut; //FilePos
  FilePos.open(FileNameOut.c_str());
  InputFile >> FileNameOut; //FileTime
  FileTime.open(FileNameOut.c_str());
  InputFile >> FileNameOut; //FileEnergy;
  FileEnergy.open(FileNameOut.c_str());
  InputFile >> FileNameOut; //FileMom;
  FileMom.open(FileNameOut.c_str());
  InputFile >> FileNameOut; //FileSor;
  FileSor.open(FileNameOut.c_str());
  InputFile >> FileNameOut; //FileVar;
  FileVar.open(FileNameOut.c_str());
  InputFile >> FileNameOut; //FilePow;
  FilePow.open(FileNameOut.c_str());
  InputFile.close();
  */
}
