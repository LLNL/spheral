#ifndef _File_Defined_
#define _File_Defined_
namespace FractalSpace
{
  class File
  {
  public:
    string BaseDirectory;
    string Directory;
    string RUN;
    ofstream FileHypre;
    ofstream FileMisc;
    ofstream FileFractalMemory;
    ofstream FileParticle;
    ofstream FileGroup;
    ofstream FileFractal;
    ofstream FilePoint;
    ofstream FileDau;
    ofstream FileEnergy;
    ofstream FileMom;
    ofstream FileForce;
    ofstream FilePow;
    ofstream FileVar;
    ofstream FileSor;
    ofstream FilePos;
    ofstream FileTimeLev;
    ofstream FileTime;
    File()
    {
    }
    File(const string& basedirectory,const int& Rank,const string& Run)
    {
      string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      BaseDirectory=basedirectory;
      cout << BaseDirectory << endl;
      stringstream streamRank;
      streamRank << Rank;
      string sRank=streamRank.str();
      cout << sRank << endl;
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+extras[ni]+"_"+sRank+"/";
	  cout << Directory << endl;
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  cout << cDirectory << endl;
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  cout << "test " << testit << " " << errno << endl;
	  if(testit == 0)
	    break;
	  assert(ni < 25);
	}
      generate_file(FileHypre,Directory+"hypre.d");
      generate_file(FileMisc,Directory+"misc.d");
      generate_file(FileFractalMemory,Directory+"fractal_memory.d");
      generate_file(FileParticle,Directory+"particle.d");
      generate_file(FileGroup,Directory+"group.d");
      generate_file(FileFractal,Directory+"fractal.d");
      generate_file(FilePoint,Directory+"point.d");
      generate_file(FileDau,Directory+"dau.d");
      generate_file(FileEnergy,Directory+"energy.d");
      generate_file(FileMom,Directory+"mom.d");
      generate_file(FileForce,Directory+"force.d");
      generate_file(FilePow,Directory+"pow.d");
      generate_file(FileVar,Directory+"var.d");
      generate_file(FileSor,Directory+"sor.d");
      generate_file(FilePos,Directory+"pos.d");
      generate_file(FileTimeLev,Directory+"time_lev.d");
      generate_file(FileTime,Directory+"time.d");

      FileMom.precision(5);
      FileForce.precision(5);
      FilePow.precision(4);
      FileVar.precision(4);
      FilePos.precision(7);
      FileTimeLev.precision(2);
      FileEnergy.precision(5);
      FileTime.precision(2);
      cout << "finished file " << endl;
    }
    void generate_file(ofstream& File,const string& sfile)
    {
      char cfile[100];
      size_t string_length=sfile.copy(cfile,1000,0);
      cfile[string_length]='\0';
      File.open(cfile);
      File << "now open for business " << endl;
    }
    void note(const bool& doit,const string& stringy)
    {
      if(doit)
	FileFractal << " " << stringy << endl;
      else
	cout << " " << stringy << endl;
    }
    ~File()
    {
      FileMisc.close();
      FileFractalMemory.close();
      FileParticle.close();
      FileGroup.close();
      FileFractal.close();
      FilePoint.close();
      FileDau.close();
      FileEnergy.close();
      FileMom.close();
      FileForce.close();
      FilePow.close();
      FileVar.close();
      FileSor.close();
      FilePos.close();
      FileTimeLev.close();
      FileTime.close();
    }
  };
}
#endif
