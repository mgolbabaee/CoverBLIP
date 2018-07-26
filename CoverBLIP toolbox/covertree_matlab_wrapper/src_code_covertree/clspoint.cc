#include "../src_hdr/stack.h"
//#define NDEBUG
#include<stdio.h>
#include<stdlib.h>
#include <complex>
#include<assert.h>
#include<math.h>
#include <string.h>
#include <iostream>
#include "mex.h"

using namespace std;

struct point_str {
     std::complex<float> pv;
     int label;
};
typedef point_str* point;

#define  PT  clspoint
class PT
{
    public:
        PT (); // constructor
        ~PT (); // DE-structor
        double distance(point v1, point v2, double upper_bound, int *counter);
        v_array<point > parse_points(float *input_r, float *input_i, int m, int n);
        v_array<point>  parse_points_chop(float *input_r, float *input_i, int chop_size, int m, int n);
        void print(point &p);
        void output(double *output_real, double *output_imag, point &p);
        void freeup (v_array<point>& pts);
        void freeup2 (v_array<v_array<point>>& pts);
        
    private:    
        const int batch = 120;//must be a multiple of 8
        int point_len_fixed;
        int point_len;
 };

/***************************************************************************************/
/* constructor */
PT::PT ()
{
    // init
    point_len_fixed = 0;
    point_len = 0;
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor */
PT::~PT ()
{
    // delete the classpoint
    point_len_fixed = 0;
    point_len = 0;
}
/***************************************************************************************/


double PT::distance(point p1, point p2, double upper_bound, int *counter)
{
  double sum = 0.;
  point_str *end = p1 + point_len_fixed;
  upper_bound *= upper_bound;
  for (point_str *batch_end = p1 + batch; batch_end <= end; batch_end += batch)
    {
      for (; p1 != batch_end; p1+=2, p2+=2)
	{
	  complex<float> d1 = p1->pv - p2->pv;
	  complex<float> d2 = (p1+1)->pv - (p2+1)->pv;
	  double d1_abs = pow(d1.real(),2) + pow(d1.imag(),2);
	  double d2_abs = pow(d2.real(),2) + pow(d2.imag(),2);
	  sum = sum + d1_abs + d2_abs;
	}
      if (sum > upper_bound)
	return sqrt(sum);
    }
  for (; p1 != end; p1+=2, p2+=2)
	{
	  complex<float> d1 = p1->pv - p2->pv;
	  complex<float> d2 = (p1+1)->pv - (p2+1)->pv;
	  double d1_abs = pow(d1.real(),2) + pow(d1.imag(),2);
	  double d2_abs = pow(d2.real(),2) + pow(d2.imag(),2);
	  sum = sum + d1_abs + d2_abs;
	}
  counter[0] = counter[0] + 1;
  return sqrt(sum);
}

v_array<point > PT::parse_points(float *input_r, float *input_i, int m, int n)
{
  v_array<point > parsed;
  v_array<point_str> p;
  point_str f;
  int i,j;
  if(input_i)
  {
  for(i = 0; i < m; i++)
  {
    int lab = i+1;
    for(j = 0; j < n; j++)
    {
      f.pv = complex<float>(input_r[i+j*m],input_i[i+j*m]);
      f.label = lab;
      push(p,f);
    }    
    if(i == 0)
        point_len = p.index;
     if (p.index %8 > 0)
	for (int k = 8 - p.index %8; k> 0; k--)
	 { f.pv = complex<float>((float) 0.,(float) 0.); f.label = lab; push(p,f);}
      point_str *new_p;
      posix_memalign((void **)&new_p, 16, p.index*sizeof(point_str));
      memcpy(new_p,p.elements,sizeof(point_str)*p.index);
      if (point_len_fixed > 0 && point_len_fixed != p.index)
	{
	  mexPrintf("Can't handle vectors of differing length, bailing\n");
	  system(0);
	}      

      point_len_fixed = p.index;
      p.index = 0;
      push(parsed,new_p);  
      //free(new_p);
  }
  }
  else
  {
  for(i = 0; i < m; i++)
  {
    int lab = i+1;
    for(j = 0; j < n; j++)
    {
      f.pv = complex<float>(input_r[i+j*m],(double) 0.);
      f.label = lab;
      push(p,f);
    }    
    if(i == 0)
        point_len = p.index;
     if (p.index %8 > 0)
	for (int k = 8 - p.index %8; k> 0; k--)
	 { f.pv = complex<float>((float) 0.,(float) 0.); f.label = lab; push(p,f);}
      point_str *new_p;
      posix_memalign((void **)&new_p, 16, p.index*sizeof(point_str));
      memcpy(new_p,p.elements,sizeof(point_str)*p.index);
      if (point_len_fixed > 0 && point_len_fixed != p.index)
	{
	  printf("Can't handle vectors of differing length, bailing\n");
	  system(0);
	}      

      point_len_fixed = p.index;
      p.index = 0;
      push(parsed,new_p);  
      //free(new_p);
  }
  }
  free(p.elements);
  return parsed;
}

v_array<point > PT::parse_points_chop(float *input_r, float *input_i, int chop_size, int m, int n)
{
  v_array<point > parsed;
  v_array<point_str> p;
  point_str f;
  int i,j;
  if(input_i)
  {
  for(i = 0; i < chop_size; i++)
  {
    int lab = i+1;
    for(j = 0; j < n; j++)
    {
      f.pv = complex<float>(input_r[i+j*m],input_i[i+j*m]);
      f.label = lab;
      push(p,f);
    }    
    if(i == 0)
        point_len = p.index;
     if (p.index %8 > 0)
	for (int k = 8 - p.index %8; k> 0; k--)
	 { f.pv = complex<float>((float) 0.,(float) 0.); f.label = lab; push(p,f);}
      point_str *new_p;
      posix_memalign((void **)&new_p, 16, p.index*sizeof(point_str));
      memcpy(new_p,p.elements,sizeof(point_str)*p.index);
      if (point_len_fixed > 0 && point_len_fixed != p.index)
	{
	  mexPrintf("Can't handle vectors of differing length, bailing\n");
	  system(0);
	}      

      point_len_fixed = p.index;
      p.index = 0;
      push(parsed,new_p);  
      //free(new_p);
  }
  }
  else
  {
  for(i = 0; i < chop_size; i++)
  {
    int lab = i+1;
    for(j = 0; j < n; j++)
    {
      f.pv = complex<float>(input_r[i+j*m],(double) 0.);
      f.label = lab;
      push(p,f);
    }    
    if(i == 0)
        point_len = p.index;
     if (p.index %8 > 0)
	for (int k = 8 - p.index %8; k> 0; k--)
	 { f.pv = complex<float>((float) 0.,(float) 0.); f.label = lab; push(p,f);}
      point_str *new_p;
      posix_memalign((void **)&new_p, 16, p.index*sizeof(point_str));
      memcpy(new_p,p.elements,sizeof(point_str)*p.index);
      if (point_len_fixed > 0 && point_len_fixed != p.index)
	{
	  printf("Can't handle vectors of differing length, bailing\n");
	  system(0);
	}      

      point_len_fixed = p.index;
      p.index = 0;
      push(parsed,new_p);  
      //free(new_p);
  }
  }
  free(p.elements);
  return parsed;
}

void PT::print(point &p)
{
  for (int i = 0; i<point_len; i++)
      mexPrintf("%f + j%f ",p[i].pv.real(), p[i].pv.imag());
  mexPrintf("\n");
}

void PT::output(double *output_real, double *output_imag, point &p)
{
//  printf("point_len = %d", point_len);
  for (int i = 0; i<point_len; i++)
  {
      output_real[i] = p[i].pv.real();
      output_imag[i] = p[i].pv.imag();
  }
}

void PT::freeup (v_array<point>& pts)
{
    for(int i = 0; i < pts.index; i++)
        if(pts[i])
            free(pts[i]);
}

void PT::freeup2 (v_array<v_array<point>>& pts)
{
    for(int i = 0; i < pts.index; i++)
        free(pts[i].elements);
    free(pts.elements);
}

