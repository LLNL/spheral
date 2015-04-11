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
    ofstream DUMPS;
    ofstream FileFractal;
    ofstream FileEnergy;
    FILE* PFHypre;
    FILE* PFFractalMemory;
    FILE* PFDau;
    FILE* PFTimeLev;
    FILE* PFTime;
    FILE* PFPos;
    FILE* PFSurface;
    File()
    {}
    File(const string& basedirectory,const int& Rank,const string& Run)
    {
      string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      BaseDirectory=basedirectory;
      if(Rank == 0)
	cerr << BaseDirectory << "\n";
      stringstream streamRank;
      streamRank << Rank;
      string sRank=streamRank.str();
      if(Rank == 0)
	cerr << sRank << "\n";
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+extras[ni]+"_"+sRank+"/";
	  if(Rank == 0)
	    cerr << Directory << "\n";
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  if(Rank == 0)
	    cerr << cDirectory << "\n";
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  if(Rank == 0)
	    cerr << "test " << testit << " " << errno << "\n";
	  if(testit == 0)
	    break;
	  assert(ni < 25);
	}
      generate_file(DUMPS,Directory+"DUMPS.d");
      generate_file(FileFractal,Directory+"fractal.d");
      generate_file(FileEnergy,Directory+"energy.d");
      PFFractalMemory=generate_FILE(Directory+"fractal_memory.d",100000);
      PFDau=PFFractalMemory;
      PFHypre=PFFractalMemory;
      PFTime=generate_FILE(Directory+"time.d",10000);
      PFTimeLev=PFTime;
      PFPos=generate_FILE(Directory+"pos.d",1000000);
      PFSurface=generate_FILE(Directory+"surface.d",1000000);
      DUMPS.precision(5);
      FileEnergy.precision(5);
      if(Rank == 0)
	cerr << "finished file " << "\n";
    }
    File(string& basedirectory,const int& FractalNodes,const int& Rank,const string& Run)
    {
      string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      //      cerr << " really basedirectory a " << basedirectory << "\n";
      stringstream ssFN;
      ssFN << FractalNodes;
      string stringFN=ssFN.str();
      basedirectory.append(stringFN);
      basedirectory.append("/");
      if(Rank == 0)
	cerr << " really basedirectory b " << basedirectory << "\n";


      BaseDirectory=basedirectory;
      //      cerr << " new BaseDirectory " << BaseDirectory << "\n";
      stringstream streamRank;
      streamRank << Rank;
      string sRank=streamRank.str();
      if(Rank == 0)
	cerr << sRank << "\n";
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+extras[ni]+"_"+sRank+"/";
	  //	  cerr << " Really Directory " << Directory << "\n";
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  if(Rank == 0)
	    cerr << cDirectory << "\n";
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  if(Rank == 0)
	    cerr << "test " << testit << " " << errno << "\n";
	  if(testit == 0)
	    break;
	  assert(ni < 25);
	}
      generate_file(DUMPS,Directory+"DUMPS.d");
      generate_file(FileFractal,Directory+"fractal.d");
      generate_file(FileEnergy,Directory+"energy.d");
      PFFractalMemory=generate_FILE(Directory+"fractal_memory.d",100000);
      PFDau=PFFractalMemory;
      PFHypre=PFFractalMemory;
      PFTime=generate_FILE(Directory+"time.d",10000);
      PFTimeLev=PFTime;
      PFPos=generate_FILE(Directory+"pos.d",1000000);
      PFSurface=generate_FILE(Directory+"surface.d",1000000);
      DUMPS.precision(5);
      FileEnergy.precision(5);
      if(Rank == 0)
	cerr << "finished file " << "\n";
    }
    ~File()
    {
      DUMPS.close();
      FileFractal.close();
      FileEnergy.close();
      fclose(PFFractalMemory);
      fclose(PFTime);
      fclose(PFPos);
      fclose(PFSurface);
    }
    void generate_file(ofstream& File,const string& sfile);
    FILE* generate_FILE(const string& sfile,int BSIZE);
    void generate_file_in(ifstream& File,const string& sfile);
    void note(const bool& doit,const string& stringy);
    void FlushAll();
  };
}
#endif
