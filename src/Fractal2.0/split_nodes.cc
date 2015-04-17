#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void split_nodes(int FR,int& FR0,int& FR1,int& FR2)
  {
    bool easy=false;
    vector <int>divs;
    factors(FR,divs,easy);
    if(easy)
      {
	FR0=divs[0];
	FR1=FR0;
	FR2=FR0;
	return;
      }
    int sumF=10*FR;
    int numfactors=divs.size();
    vector <int>Nodes(3);
    vector <int>threepow(numfactors);
    threepow[0]=1;
    for(int num=1;num<numfactors;num++)
      threepow[num]=threepow[num-1]*3;
    int options=threepow[numfactors-1]*3;
    for(int nopt=0;nopt<options;nopt++)
      {
	Nodes.assign(3,1);
	for(int num=0;num<numfactors;num++)
	  {
	    int which=(nopt/threepow[num]) % 3;
	    Nodes[which]*=divs[num];
	  }
	int sum3=Nodes[0]+Nodes[1]+Nodes[2];
	if(sum3 < sumF)
	  {
	    sumF=sum3;
	    FR0=Nodes[0];
	    FR1=Nodes[1];
	    FR2=Nodes[2];
	  }
      }
  }
  void factors(int FR,vector <int>& divs,bool& easy)
  {
    //  cerr << " enter factors " << FR << endl;
    divs.clear();
    int third=pow(static_cast<double>(FR)+0.5,1.0/3.0);
    if(third*third*third == FR)
      {
	divs.push_back(third);
	easy=true;
	return;
      }
    divs.push_back(1);
    int diva=FR;
    int divisor=2;
    while(divisor <= FR && diva > 1)
      {
	if(diva % divisor == 0)
	  {
	    divs.push_back(divisor);
	    diva/=divisor;
	  }
	else
	  divisor++;
	//      cerr << " In factors " << FR << " " << diva << " " << divisor << " " << endl;
      }
  }
}
