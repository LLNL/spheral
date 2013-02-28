#ifndef _Misc_Defined_
#define _Misc_Defined_
namespace FractalSpace
{
  class Misc
  {
    bool debug;
  public:
    Group* p_group_0;
    int zoom;
    int grid_multiply;
    static int dim0;
    static int dim1;
    static int dim2;
    Misc()
    {
      assert(this);
      debug=false;
      cout << "Making Misc " << this << endl;
    }
    ~Misc(){
      cout << "Ending Misc " << this << endl;
    }
    bool get_debug()
    {
      return debug;
    }
    void set_debug(bool& d)
    {
      debug=d;
    }
    static int nr(const int& i,const int& j,const int& k, const int&m)
    {
      return i+(j+k*m)*m;
    }
    static int pow(const int&x,const int& y)
    {
      assert(y >= 0);
      int i=1;
      for(int ii=0;ii < y;++ii)
	i*=x;
      return i;
    }
    template <class T>  static T pow2(const T& x)
    {
      return x*x;
    }
    template <class T>  static T pow3(const T& x)
    {
      return x*x*x;
    }
    template <class T> static void add_dens(vector <T>& dens,const T& dm, const T& d_x, 
					    const T& d_y, const T& d_z)
    {
      assert(abs(d_x-0.5) <= 0.5);
      assert(abs(d_y-0.5) <= 0.5);
      assert(abs(d_z-0.5) <= 0.5);
      T d_1=(1.0-d_x)*(1.0-d_y);
      T d_2=d_x*(1.0-d_y);
      T d_3=(1.0-d_x)*d_y;
      T d_4=d_x*d_y;
      T d_z_1_dm=(1.0-d_z)*dm;
      T d_z_dm=d_z*dm;
      //
      dens[0]+=d_1*d_z_1_dm;
      dens[1]+=d_2*d_z_1_dm;
      dens[2]+=d_3*d_z_1_dm;
      dens[3]+=d_4*d_z_1_dm;
      dens[4]+=d_1*d_z_dm;
      dens[5]+=d_2*d_z_dm;
      dens[6]+=d_3*d_z_dm;
      dens[7]+=d_4*d_z_dm;
    }
    template <class T> static T sum_prod(const int& n1,const int& n2,const int& n3, vector <T>& x, vector <T>& y)
    {
      T  sum=0.0;
      for (int n=n1 ; n <= n2; n+=n3)
	{
	  sum+=x[n]*y[n];
	}
      return sum;
    }
    template <class T>   static void sum_prod(const int& n1,const int& n2,const int& n3, vector <T>& sum_4,vector <T>& x, 
					      vector <T>& a, vector <T>& b, vector <T>& c, vector <T>& d)
    {
      sum_4[0]=0.0;
      sum_4[1]=0.0;
      sum_4[2]=0.0;
      sum_4[3]=0.0;
      for (int n=n1 ; n <= n2; n+=n3)
	{
	  sum_4[0]+=x[n]*a[n];
	  sum_4[1]+=x[n]*b[n];
	  sum_4[2]+=x[n]*c[n];
	  sum_4[3]+=x[n]*d[n];
	}
    }
    template <class T>   static void sum_prod(const int& n1,const int& n2,const int& n3, vector <T>& sum_4,vector <T>& x, 
					      vector <T>& a, vector <T>& b, vector <T>& c, vector <T>& d, vector <T>& e, vector <T>& f)
    {
      sum_4[0]=0.0;
      sum_4[1]=0.0;
      sum_4[2]=0.0;
      sum_4[3]=0.0;
      sum_4[4]=0.0;
      sum_4[5]=0.0;
      for (int n=n1 ; n <= n2; n+=n3)
	{
	  sum_4[0]+=x[n]*a[n];
	  sum_4[1]+=x[n]*b[n];
	  sum_4[2]+=x[n]*c[n];
	  sum_4[3]+=x[n]*d[n];
	  sum_4[4]+=x[n]*e[n];
	  sum_4[5]+=x[n]*f[n];
	}
    }
    template <class T>  static void sum_prod_p_sharp(const int& n1,const int& n2,const int& n3, vector <T>& sum_4,
						     vector <T>& w_p,vector <T>& w_x,vector <T>& w_y,vector <T>& w_z, 
						     vector <T>& a)
    {
      sum_4[0]=0.0;
      sum_4[1]=0.0;
      sum_4[2]=0.0;
      sum_4[3]=0.0;
      for (int n=n1 ; n <= n2; n+=n3)
	{
	  sum_4[0]+=w_p[n]*a[n];
	  sum_4[1]+=w_x[n]*a[n];
	  sum_4[2]+=w_y[n]*a[n];
	  sum_4[3]+=w_z[n]*a[n];
	}
    }
    template <class T> static void set_weights(vector <T>& weights,const T& d_x,const T& d_y,const T& d_z)
    {
      T d_z_1=1.0-d_z;
      weights[0]=(1.0-d_x)*(1.0-d_y);
      weights[1]=d_x*(1.0-d_y);
      weights[2]=(1.0-d_x)*d_y;
      weights[3]=d_x*d_y;
      weights[4]=weights[0]*d_z;
      weights[5]=weights[1]*d_z;
      weights[6]=weights[2]*d_z;
      weights[7]=weights[3]*d_z;
      weights[0]*=d_z_1;
      weights[1]*=d_z_1;
      weights[2]*=d_z_1;
      weights[3]*=d_z_1;
    }
    template <class T> static void set_weights(vector <T>& weights_p,vector <T>& weights_x,
					       vector <T>& weights_y,vector <T>& weights_z,
					       const T& d_x,const T& d_y,const T& d_z)
    {
      T d_z_1=1.0-d_z;
      weights_p[0]=(1.0-d_x)*(1.0-d_y);
      weights_p[1]=d_x*(1.0-d_y);
      weights_p[2]=(1.0-d_x)*d_y;
      weights_p[3]=d_x*d_y;
      weights_p[4]=weights_p[0]*d_z;
      weights_p[5]=weights_p[1]*d_z;
      weights_p[6]=weights_p[2]*d_z;
      weights_p[7]=weights_p[3]*d_z;
      weights_p[0]*=d_z_1;
      weights_p[1]*=d_z_1;
      weights_p[2]*=d_z_1;
      weights_p[3]*=d_z_1;

      weights_x[0]=(weights_p[0]+weights_p[1]);
      weights_x[1]=-weights_x[0];
      weights_x[2]=(weights_p[2]+weights_p[3]);
      weights_x[3]=-weights_x[2];
      weights_x[4]=(weights_p[4]+weights_p[5]);
      weights_x[5]=-weights_x[4];
      weights_x[6]=(weights_p[6]+weights_p[7]);
      weights_x[7]=-weights_x[6];

      weights_y[0]=(weights_p[0]+weights_p[2]);
      weights_y[2]=-weights_y[0];
      weights_y[1]=(weights_p[1]+weights_p[3]);
      weights_y[3]=-weights_y[1];
      weights_y[4]=(weights_p[4]+weights_p[6]);
      weights_y[6]=-weights_y[4];
      weights_y[5]=(weights_p[5]+weights_p[7]);
      weights_y[7]=-weights_y[5];

      weights_z[0]=(weights_p[0]+weights_p[4]);
      weights_z[4]=-weights_z[0];
      weights_z[1]=(weights_p[1]+weights_p[5]);
      weights_z[5]=-weights_z[1];
      weights_z[2]=(weights_p[2]+weights_p[6]);
      weights_z[6]=-weights_z[2];
      weights_z[3]=(weights_p[3]+weights_p[7]);
      weights_z[7]=-weights_z[3];
    }
    template <class T>  static T sinc_2(const T& x)
    {
      if(abs(x) > 0.0001)
	{
	  T s=sin(x)/x;
	  return s*s;
	}
      return 1.0-x*x/3.0;
    } 
    template <class T>  static T square_filter(const T& x)
    {
      if(abs(x) > 0.0001)
	{
	  T s=3.0*(sin(x)-x*cos(x))/(x*x*x);
	  return s*s;
	}
      return 1.0-x*x*0.2;
    }
    template <class T> static void my_assign(vector <T>& vector1,const int& itr1_begin ,const int& itr1_d,const T& value)
    {
      for(int itr1=itr1_begin;itr1<itr1_begin+itr1_d;itr1++)
	{
	  vector1[itr1]=value;
	}
    }
    template <class T> void per_box(vector <T>& box,const T& length)
    {
      unsigned int bs=box.size();
      if(bs == 6)
	{
	  for(int ni2=0;ni2<6;ni2+=2)
	    {
	      db=box[ni2+1]-box[ni2];
	      box[ni2]=(box[ni2]+length) % length;
	      box[ni2+1]=box[ni2]+db;
	    }
	}
      else if(bs == 3)
	{
	  box[0]=(box[0]+length) % length;
	  box[1]=(box[1]+length) % length;
	  box[2]=(box[2]+length) % length;
	}
      else
	assert(0);
    }
    template <class T> void per_box(vector <T>& posa,vector <T>& posb,const T& length)
    {
      unsigned int ps=posa.size();
      assert(ps == posb.size());
      T dp;
      for(unsigned int ni=0;ni<ps;ni++)
	{
	  dp=posb[ni]-posa[ni];
	  posa[ni]=(posa[ni]+length) % length;
	  posb[ni]=posa[ni]+dp;
	}
    }
    template <class T> static void copy_vector(int& itr1,vector <T>& vector1,const vector <T>& vector2,const int& itr2_begin ,const int& itr2_d)
    {
      for(int i=itr2_begin;i<itr2_begin+itr2_d;i++)
	{
	  vector1[itr1]=vector2[i];
	  itr1++;
	}
    }
    template <class T> static void copy_vector_non_zero(int& itr1,vector <T>& vector1,const vector <T>& vector2,const int& itr2_begin ,const int& itr2_d)
    {
      for(int itr2=itr2_begin;itr2<itr2_begin+itr2_d;itr2++)
	{
	  if(vector2[itr2] != 0)
	    {
	      vector1[itr1]=vector2[itr2];
	      itr1++;
	    }
	}
    }
    template <class T> static void vector_print(const vector <T>& vec,ofstream& FILE)
    {
      int j=vec.size();
      for(int i=0;i<j;i++)
	FILE << vec[i] << " " ;
      FILE << endl;
    }
    template <class T> static void vector_print(const vector <T>& veca,const vector <T>& vecb,ofstream& FILE)
    {
      int j=veca.size();
      for(int i=0;i<j;i++)
	FILE << veca[i] << " " ;
      int j=vecb.size();
      for(int i=0;i<j;i++)
	FILE << vecb[i] << " " ;
      FILE << endl;
    }
    template <class T> static void vector_print(const vector <T>& veca,const vector <T>& vecb,const vector <T>& vecc,ofstream& FILE)
    {
      int j=veca.size();
      for(int i=0;i<j;i++)
	FILE << veca[i] << " " ;
      int j=vecb.size();
      for(int i=0;i<j;i++)
	FILE << vecb[i] << " " ;
      int j=vecc.size();
      for(int i=0;i<j;i++)
	FILE << vecc[i] << " " ;
      FILE << endl;
    }
  };
}
#endif
