#ifndef _Misc_Defined_
#define _Misc_Defined_

#include <vector>
#include <fstream>

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
    }
    ~Misc()
    {
    }
    bool get_debug() const;
    void set_debug(bool& d);
    static int coordinate(std::vector <int>& pos,std::vector <int>& Box,int spacing)
    {
      int nx=pos[0]-Box[0];
      int ny=pos[1]-Box[2];
      int nz=pos[2]-Box[4];
      int nxt=(Box[1]-Box[0]+1)/spacing;
      int nyt=(Box[3]-Box[2]+1)/spacing;
      return (nx+nxt*(ny+nz*nyt))/spacing;
    }
     template <class T> void plus(const std::vector <T>& vin,T addit,std::vector <T>& vout)
    {
      vout=vin;
      plus(vout,addit);
    }
    template <class T> static void plus(std::vector <T>& vect,T add)
    {
      for(auto &v : vect)
	v+=add;
    }
    template <class T> static void times(const std::vector <T>& vin,const T mult,std::vector <T>& vout)
    {
      vout=vin;
      times(vout,mult);
    }
    template <class T> static void times(std::vector <T>& vect,const T mult)
    {
      for(auto &v : vect)
	v*=mult;
    }
    template <class T> static void divide(const std::vector <T>& vin,const T divisor,std::vector <T>& vout)
    {
      vout=vin;
      divide(vout,divisor);
    }
    template <class T> static void divide(std::vector <T>& vect,const T divisor)
    {
      for(auto &v : vect)
	v/=divisor;
    }
    template <class T> static T nr(const T& i,const T& j,const T& k, const T&m)
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
    template <class T> static T pow2(const T& x)
    {
      return x*x;
    }
    template <class T> static T pow3(const T& x)
    {
      return x*x*x;
    }
    template <class T> static void add_dens(std::vector <T>& dens,const T& dm, T& d_x, 
					    T& d_y, T& d_z)
    {
      //       assert(abs(d_x-0.5) <= 0.5);
      //       assert(abs(d_y-0.5) <= 0.5);
      //       assert(abs(d_z-0.5) <= 0.5);
      if(abs(d_x-0.5) >= 0.5)
	{
	  std::cout << "dx error " << abs(d_x-0.5)-0.5 << "\n";
	  if(d_x > 1.0)
	    d_x=0.9999999;
	  else
	    d_x=1.0e-6;
	}
      if(abs(d_y-0.5) >= 0.5)
	{
	  std::cout << "dy error " << abs(d_y-0.5)-0.5 << "\n";
	  if(d_y > 1.0)
	    d_y=0.9999999;
	  else
	    d_y=1.0e-6;
	}
      if(abs(d_z-0.5) >= 0.5)
	{
	  std::cout << "dz error " << abs(d_z-0.5)-0.5 << "\n";
	  if(d_z > 1.0)
	    d_z=0.9999999;
	  else
	    d_z=1.0e-6;
	}
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
    template <class T> static T sum_prod(const int& n1,const int& n2,const int& n3, std::vector <T>& x, std::vector <T>& y)
    {
      T  sum=0.0;
      for (int n=n1 ; n <= n2; n+=n3)
	{
	  sum+=x[n]*y[n];
	}
      return sum;
    }
    template <class T> static void sum_prod(const int& n1,const int& n2,const int& n3, std::vector <T>& sum_4,std::vector <T>& x, 
					    std::vector <T>& a, std::vector <T>& b, std::vector <T>& c, std::vector <T>& d)
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
    template <class T> static void sum_prod(const int& n1,const int& n2,const int& n3, std::vector <T>& sum_4,std::vector <T>& x, 
					    std::vector <T>& a, std::vector <T>& b, std::vector <T>& c, std::vector <T>& d, std::vector <T>& e, std::vector <T>& f)
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
    template <class T> static void sum_prod_p_sharp(const int& n1,const int& n2,const int& n3, std::vector <T>& sum_4,
						    std::vector <T>& w_p,std::vector <T>& w_x,std::vector <T>& w_y,std::vector <T>& w_z, 
						    std::vector <T>& a)
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
    template <class T> static void set_weights(std::vector <T>& weights,const T& d_x,const T& d_y,const T& d_z)
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
    template <class T> static void set_weights(std::vector <T>& weights_p,std::vector <T>& weights_x,
					       std::vector <T>& weights_y,std::vector <T>& weights_z,
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
    template <class T> static T sinc_2(const T& x)
    {
      if(abs(x) > 0.0001)
	{
	  T s=sin(x)/x;
	  return s*s;
	}
      return 1.0-x*x/3.0;
    } 
    template <class T> static T square_filter(const T& x)
    {
      if(abs(x) > 0.0001)
	{
	  T s=3.0*(sin(x)-x*cos(x))/(x*x*x);
	  return s*s;
	}
      return 1.0-x*x*0.2;
    }
    template <class T> static void my_assign(std::vector <T>& vector1,const int& itr1_begin ,const int& itr1_d,const T& value)
    {
      for(int itr1=itr1_begin;itr1<itr1_begin+itr1_d;itr1++)
	{
	  vector1[itr1]=value;
	}
    }
    template <class T> void per_box(std::vector <T>& box,const T& length)
    {
      unsigned int bs=box.size();
      if(bs == 6)
	{
	  for(int ni2=0;ni2<6;ni2+=2)
	    {
	      int db=box[ni2+1]-box[ni2];
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
    template <class T> void per_box(std::vector <T>& posa,std::vector <T>& posb,const T& length)
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
    template <class T> static void copy_vector(int& itr1,std::vector <T>& vector1,const std::vector <T>& vector2,const int& itr2_begin ,const int& itr2_d)
    {
      for(int i=itr2_begin;i<itr2_begin+itr2_d;i++)
	{
	  vector1[itr1]=vector2[i];
	  itr1++;
	}
    }
    template <class T> static void copy_vector_non_zero(int& itr1,std::vector <T>& vector1,const std::vector <T>& vector2,const int& itr2_begin ,const int& itr2_d)
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
    template <class T> static void vector_print(const std::vector <T>& vec,std::ofstream& FILE)
    {
      int j=vec.size();
      for(int i=0;i<j;i++)
	FILE << vec[i] << " " ;
      FILE << "\n";
    }
    template <class T> static void vector_print(const std::vector <T>& veca,const std::vector <T>& vecb,std::ofstream& FILE)
    {
      int j=veca.size();
      for(int i=0;i<j;i++)
	FILE << veca[i] << " " ;
      j=vecb.size();
      for(int i=0;i<j;i++)
	FILE << vecb[i] << " " ;
      FILE << "\n";
    }
    template <class T> static void vector_print(const std::vector <T>& veca,const std::vector <T>& vecb,const std::vector <T>& vecc,std::ofstream& FILE)
    {
      int j=veca.size();
      for(int i=0;i<j;i++)
	FILE << veca[i] << " " ;
      j=vecb.size();
      for(int i=0;i<j;i++)
	FILE << vecb[i] << " " ;
      j=vecc.size();
      for(int i=0;i<j;i++)
	FILE << vecc[i] << " " ;
      FILE << "\n";
    }
    template <class T> static void sum_up(T& sum,std::vector <T>& values,int first,int last,const int stride)
    {
      sum=0;
      if((last-first)*stride <= 0)
	return;
      while(first < last)
	{
	  sum+=values[first];
	  first+=stride;
	}
      return;
    }
    // template <class T> void zero_shrink_std::vector(std::vector <T>& vec,int size)
    // {
    //   vec.clear();
    //   vec.resize(size);
    //   vec.shrink_to_fit();
    // }
    // template <class T> void shrink_std::vectors(std::vector <std::vector <T> >& vec)
    // {
    //   vec.shrink_to_fit();
    //   for(std::vector <T> v : vec)
    // 	v.shrink_to_fit();
    // }
    template <class T> struct count_up
    {
      bool operator()(const T A,const T B) const
      {
	return A < B;
      }
    };
    template <class T> struct count_down
    {
      bool operator()(const T A,const T B) const
      {
	return A > B;
      }
    };
  };
}
#endif
