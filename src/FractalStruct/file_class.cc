#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void File::generate_file(ofstream& File,const string& sfile)
  {
    char cfile[100];
    size_t string_length=sfile.copy(cfile,1000,0);
    cfile[string_length]='\0';
    File.open(cfile);
    //      cerr << " newFile " << cfile << " " << &File << "\n";
    //      File << "now open for business " << "\n";
  }
  FILE* File::generate_FILE(const string& sfile,int BSIZE)
  {
    char cfile[100];
    size_t string_length=sfile.copy(cfile,1000,0);
    cfile[string_length]='\0';
    FILE* PFile=fopen(cfile,"w");
    setvbuf(PFile,NULL,_IOFBF,BSIZE);
    //      cerr << " newFILE " << cfile << " " << BSIZE << " " << PFile << "\n";
    return PFile;
  }
  void File::generate_file_in(ifstream& File,const string& sfile)
  {
    char cfile[100];
    size_t string_length=sfile.copy(cfile,1000,0);
    cfile[string_length]='\0';
    File.open(cfile);
    //      File << "now open for business " << "\n";
  }
  void File::note(const bool& doit,const string& stringy)
  {
    if(doit)
      FileFractal << " " << stringy << "\n";
    else
      cerr << " " << stringy << "\n";
  }
  void File::FlushAll()
  {
    FileFractal.flush();
    FileEnergy.flush();
    cerr.flush();
    fflush(PFFractalMemory);
    fflush(PFTime);
    fflush(PFPos);
    fflush(PFBox);
    fflush(PFSurface);
  }
}
