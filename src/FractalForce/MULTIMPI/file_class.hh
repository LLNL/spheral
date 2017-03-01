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
    //    ofstream FileFFT;
    //    ofstream FileHypreTime;
    //    ofstream FileMisc;
    //    ofstream FileParticle;
    //    ofstream FileGroup;
    ofstream FileFractal;
    //    ofstream FilePoint;
    ofstream FileEnergy;
    //    ofstream FileMom;
    //    ofstream FileForce;
    //    ofstream FilePow;
    //    ofstream FileVar;
    //    ofstream FileSor;
    FILE* PFHypre;
    FILE* PFFractalMemory;
    FILE* PFDau;
    FILE* PFTimeLev;
    FILE* PFTime;
    FILE* PFPos;
    //    FILE* PFTree;
    //    FILE* PFSurface;
    File()
    {
    }
    File(const string& basedirectory,const int& Rank,const string& Run)
    {
      string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      BaseDirectory=basedirectory;
      if(Rank == 0)
	cout << BaseDirectory << "\n";
      stringstream streamRank;
      streamRank << Rank;
      string sRank=streamRank.str();
      if(Rank == 0)
	cout << sRank << "\n";
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+extras[ni]+"_"+sRank+"/";
	  if(Rank == 0)
	    cout << Directory << "\n";
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  if(Rank == 0)
	    cout << cDirectory << "\n";
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  if(Rank == 0)
	    cout << "test " << testit << " " << errno << "\n";
	  if(testit == 0)
	    break;
	  assert(ni < 25);
	}
      generate_file(DUMPS,Directory+"DUMPS.d");
      //      generate_file(FileFFT,Directory+"fft.d");
      //      generate_file(FileHypreTime,Directory+"hypretime.d");
      //      generate_file(FileMisc,Directory+"misc.d");
      //      generate_file(FileParticle,Directory+"particle.d");
      //      generate_file(FileGroup,Directory+"group.d");
      generate_file(FileFractal,Directory+"fractal.d");
      //      generate_file(FilePoint,Directory+"point.d");
      generate_file(FileEnergy,Directory+"energy.d");
      //      generate_file(FileMom,Directory+"mom.d");
      //      generate_file(FileForce,Directory+"force.d");
      //      generate_file(FilePow,Directory+"pow.d");
      //      generate_file(FileVar,Directory+"var.d");
      //      generate_file(FileSor,Directory+"sor.d");
      PFFractalMemory=generate_FILE(Directory+"fractal_memory.d",100000);
      PFDau=PFFractalMemory;
      //      PFDau=generate_FILE(Directory+"dau.d",100000);
      PFHypre=PFFractalMemory;
      //      PFHypre=generate_FILE(Directory+"hypre.d",100000);
      PFTime=generate_FILE(Directory+"time.d",10000);
      PFTimeLev=PFTime;
      //      PFTimeLev=generate_FILE(Directory+"time_lev.d",10000);
      PFPos=generate_FILE(Directory+"pos.d",1000000);
      //      PFSurface=generate_FILE(Directory+"surface.d",1000000);
      //      PFTree=PFSurface;
      //      PFTree=generate_FILE(Directory+"tree.d",1000000);
      DUMPS.precision(5);
      //      FileMom.precision(5);
      //      FileForce.precision(5);
      //      FilePow.precision(4);
      //      FileVar.precision(4);
      FileEnergy.precision(5);
      if(Rank == 0)
	cout << "finished file " << "\n";
    }
    File(string& basedirectory,const int& FractalNodes,const int& Rank,const string& Run)
    {
      string extras("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      RUN=Run;
      //      cout << " really basedirectory a " << basedirectory << "\n";
      stringstream ssFN;
      ssFN << FractalNodes;
      string stringFN=ssFN.str();
      basedirectory.append(stringFN);
      basedirectory.append("/");
      if(Rank == 0)
	cout << " really basedirectory b " << basedirectory << "\n";


      BaseDirectory=basedirectory;
      //      cout << " new BaseDirectory " << BaseDirectory << "\n";
      stringstream streamRank;
      streamRank << Rank;
      string sRank=streamRank.str();
      if(Rank == 0)
	cout << sRank << "\n";
      for(int ni=0;ni<26;ni++)
	{
	  Directory=BaseDirectory+RUN+extras[ni]+"_"+sRank+"/";
	  //	  cout << " Really Directory " << Directory << "\n";
	  char cDirectory[200];
	  size_t dir_length=Directory.copy(cDirectory,1000,0);
	  cDirectory[dir_length]='\0';
	  if(Rank == 0)
	    cout << cDirectory << "\n";
	  int testit=mkdir(cDirectory,S_IRWXU|S_IWGRP|S_IXGRP);
	  if(Rank == 0)
	    cout << "test " << testit << " " << errno << "\n";
	  if(testit == 0)
	    break;
	  assert(ni < 25);
	}
      generate_file(DUMPS,Directory+"DUMPS.d");
      //      generate_file(FileFFT,Directory+"fft.d");
      //      generate_file(FileHypreTime,Directory+"hypretime.d");
      //      generate_file(FileMisc,Directory+"misc.d");
      //      generate_file(FileParticle,Directory+"particle.d");
      //      generate_file(FileGroup,Directory+"group.d");
      generate_file(FileFractal,Directory+"fractal.d");
      //      generate_file(FilePoint,Directory+"point.d");
      generate_file(FileEnergy,Directory+"energy.d");
      //      generate_file(FileMom,Directory+"mom.d");
      //      generate_file(FileForce,Directory+"force.d");
      //      generate_file(FilePow,Directory+"pow.d");
      //      generate_file(FileVar,Directory+"var.d");
      //      generate_file(FileSor,Directory+"sor.d");
      PFFractalMemory=generate_FILE(Directory+"fractal_memory.d",100000);
      PFDau=PFFractalMemory;
      //      PFDau=generate_FILE(Directory+"dau.d",100000);
      PFHypre=PFFractalMemory;
      //      PFHypre=generate_FILE(Directory+"hypre.d",100000);
      PFTime=generate_FILE(Directory+"time.d",10000);
      PFTimeLev=PFTime;
      //      PFTimeLev=generate_FILE(Directory+"time_lev.d",10000);
      PFPos=generate_FILE(Directory+"pos.d",1000000);
      //      PFSurface=generate_FILE(Directory+"surface.d",1000000);
      //      PFTree=PFSurface;
      //      PFTree=generate_FILE(Directory+"tree.d",1000000);
      DUMPS.precision(5);
      //      FileMom.precision(5);
      //      FileForce.precision(5);
      //      FilePow.precision(4);
      //      FileVar.precision(4);
      FileEnergy.precision(5);
      if(Rank == 0)
	cout << "finished file " << "\n";
    }
    void generate_file(ofstream& File,const string& sfile)
    {
      char cfile[100];
      size_t string_length=sfile.copy(cfile,1000,0);
      cfile[string_length]='\0';
      File.open(cfile);
      //      cout << " newFile " << cfile << " " << &File << "\n";
      //      File << "now open for business " << "\n";
    }
    FILE* generate_FILE(const string& sfile,int BSIZE)
    {
      char cfile[100];
      size_t string_length=sfile.copy(cfile,1000,0);
      cfile[string_length]='\0';
      FILE* PFile=fopen(cfile,"w");
      setvbuf(PFile,NULL,_IOFBF,BSIZE);
      //      cout << " newFILE " << cfile << " " << BSIZE << " " << PFile << "\n";
      return PFile;
    }
    void generate_file_in(ifstream& File,const string& sfile)
    {
      char cfile[100];
      size_t string_length=sfile.copy(cfile,1000,0);
      cfile[string_length]='\0';
      File.open(cfile);
      //      File << "now open for business " << "\n";
    }
    void note(const bool& doit,const string& stringy)
    {
      if(doit)
	FileFractal << " " << stringy << "\n";
      else
	cout << " " << stringy << "\n";
    }
    void FlushAll()
    {
      //      FileFFT.flush();
      //      FileHypreTime.flush();
      //      FileMisc.flush();
      //      FileParticle.flush();
      //      FileGroup.flush();
      FileFractal.flush();
      //      FilePoint.flush();
      FileEnergy.flush();
      //      FileMom.flush();
      //      FileForce.flush();
      //      FilePow.flush();
      //      FileVar.flush();
      //      FileSor.flush();
      cout.flush();
      //      fflush(PFHypre);
      fflush(PFFractalMemory);
      //      fflush(PFDau);
      //      fflush(PFTimeLev);
      fflush(PFTime);
      fflush(PFPos);
      //      fflush(PFSurface);
      //      fflush(PFTree);
    }
    ~File()
    {
      DUMPS.close();
      //      FileFFT.close();
      //      FileHypreTime.close();
      //      FileMisc.close();
      //      FileParticle.close();
      //      FileGroup.close();
      FileFractal.close();
      //      FilePoint.close();
      FileEnergy.close();
      //      FileMom.close();
      //      FileForce.close();
      //      FilePow.close();
      //      FileVar.close();
      //      FileSor.close();
      //      fclose(PFDau);
      //      fclose(PFHypre);
      fclose(PFFractalMemory);
      //      fclose(PFTimeLev);
      fclose(PFTime);
      fclose(PFPos);
      //      fclose(PFTree);
      //      fclose(PFSurface);
    }
  };
}
#endif
