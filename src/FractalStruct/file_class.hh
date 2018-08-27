#ifndef _File_Defined_
#define _File_Defined_

#include <iostream>
#include <fstream>
#include <sys/stat.h>

namespace FractalSpace
{
  class File
  {
  public:
    std::string BaseDirectory;
    std::string Directory;
    std::string RUN;
    std::ofstream DUMPS;
    std::ofstream FileFractal;
    std::ofstream FileEnergy;
    FILE* PFHypre;
    FILE* PFFractalMemory;
    FILE* PFDau;
    FILE* PFTimeLev;
    FILE* PFTime;
    FILE* PFPos;
    FILE* PFBox;
    FILE* PFSurface;
    File()
    {
      DUMPS.open("/dev/null");
      FileFractal.open("/dev/null");
      FileEnergy.open("/dev/null");
      PFHypre=fopen("/dev/null","w");
      PFFractalMemory=fopen("/dev/null","w");
      PFDau=fopen("/dev/null","w");
      PFTimeLev=fopen("/dev/null","w");
      PFTime=fopen("/dev/null","w");
      PFPos=fopen("/dev/null","w");
      PFBox=fopen("/dev/null","w");
      PFSurface=fopen("/dev/null","w");
    }
    File(const std::string& basedirectory,const int& Rank,const std::string& Run)
    {
      std::string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      BaseDirectory=basedirectory;
      if(Rank == 0)
        std::cerr << BaseDirectory << "\n";
      std::stringstream streamRank;
      streamRank << Rank;
      std::string sRank=streamRank.str();
      if(Rank == 0)
	std::cerr << sRank << "\n";
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+"_"+extras[ni]+"_"+sRank+"/";
	  if(Rank == 0)
	    std::cerr << Directory << "\n";
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  if(Rank == 0)
	    std::cerr << cDirectory << "\n";
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  if(Rank == 0)
	    std::cerr << "test " << testit << " " << errno << "\n";
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
      PFBox=generate_FILE(Directory+"box.d",1000000);
      PFSurface=generate_FILE(Directory+"surface.d",1000000);
      DUMPS.precision(5);
      FileEnergy.precision(5);
      if(Rank == 0)
	std::cerr << "finished file " << "\n";
    }
    File(std::string& basedirectory,const int& FractalNodes,const int& Rank,const std::string& Run)
    {
      std::string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      //      std::cerr << " really basedirectory a " << basedirectory << "\n";
      std::stringstream ssFN;
      ssFN << FractalNodes;
      std::string stringFN=ssFN.str();
      basedirectory.append(stringFN);
      basedirectory.append("/");
      if(Rank == 0)
	std::cerr << " really basedirectory b " << basedirectory << "\n";


      BaseDirectory=basedirectory;
      //      std::cerr << " new BaseDirectory " << BaseDirectory << "\n";
      std::stringstream streamRank;
      streamRank << Rank;
      std::string sRank=streamRank.str();
      if(Rank == 0)
	std::cerr << sRank << "\n";
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+"_"+extras[ni]+"_"+sRank+"/";
	  //	  std::cerr << " Really Directory " << Directory << "\n";
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  if(Rank == 0)
	    std::cerr << cDirectory << "\n";
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  if(Rank == 0)
	    std::cerr << "test " << testit << " " << errno << "\n";
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
      PFBox=generate_FILE(Directory+"box.d",1000000);
      PFSurface=generate_FILE(Directory+"surface.d",1000000);
      DUMPS.precision(5);
      FileEnergy.precision(5);
      if(Rank == 0)
	std::cerr << "finished file " << "\n";
    }
    ~File()
    {
      DUMPS.close();
      FileFractal.close();
      FileEnergy.close();
      fclose(PFFractalMemory);
      fclose(PFTime);
      fclose(PFPos);
      fclose(PFBox);
      fclose(PFSurface);
    }
    void generate_file(std::ofstream& File,const std::string& sfile);
    FILE* generate_FILE(const std::string& sfile,int BSIZE);
    void generate_file_in(std::ifstream& File,const std::string& sfile);
    void note(const bool& doit,const std::string& stringy);
    void FlushAll();
  };
}
#endif
