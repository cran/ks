#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <R_ext/Arith.h>
#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* Multivariate linear binning functions 
   translated from the Fortran code of M. Wand & T.Duong in ks < 1.8.0
   adapted from 1-d massdist.c in stats package */

/* Headers */
void massdist1d(double *x1, int *n, double *a1, double *b1, int *M1,
                double *weight, double *est);

void massdist2d(double *x1, double *x2,	int *n, 
		double *a1, double *a2, double *b1, double *b2,
		int *M1, int *M2, double *weight, double *est);

void massdist3d(double *x1, double *x2,	double *x3, int *n, 
		double *a1, double *a2, double *a3, 
		double *b1, double *b2, double *b3,
		int *M1, int *M2, int *M3, double *weight, double *est);

void massdist4d(double *x1, double *x2,	double *x3, double *x4, int *n, 
		double *a1, double *a2, double *a3, double *a4, 
		double *b1, double *b2, double *b3, double *b4, 
		int *M1, int *M2, int *M3, int *M4, double *weight, double *est);

void interp1d(double *x1, int *n, 
              double *a1, double *b1, int *M1,
              double *fun, double *est);

void interp2d(double *x1, double *x2, int *n, 
	      double *a1, double *a2, double *b1, double *b2,
	      int *M1, int *M2, double *fun, double *est);

void interp3d(double *x1, double *x2, double *x3, int *n, 
	      double *a1, double *a2, double *a3, 
	      double *b1, double *b2, double *b3,
              int *M1, int *M2, int *M3, double *fun, double *est);


/* Code */

void massdist1d(double *x1, int *n, double *a1, double *b1, int *M1,
               double *weight, double *est)
{ 
  double fx1, wi, xdelta1, xpos1;   
  int i, ix1, ixmax1, ixmin1, MM1;

  MM1 = M1[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
   
  // set all est = 0 
  for (i=0; i < MM1; i++)
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      ix1 = floor(xpos1);
      fx1 = xpos1 - ix1;
      wi = weight[i];   
      
      if(ixmin1 <= ix1 && ix1 <= ixmax1) {
	est[ix1] += wi*(1-fx1);   
	est[ix1 + 1] += wi*fx1;
      }
      else if(ix1 == -1) {
	est[0] += wi*fx1; 
      }
      else if(ix1 == ixmax1 + 1) {
	est[ix1] += wi*(1-fx1);  
      }
    }
  } 
}

void massdist2d(double *x1, double *x2,	int *n, 
		double *a1, double *a2, double *b1, double *b2,
		int *M1, int *M2, double *weight, double *est)
{
  double fx1, fx2, wi, xdelta1, xdelta2, xpos1, xpos2;   
  int i, ix1, ix2, ixmax1, ixmin1, ixmax2, ixmin2, MM1, MM2;
  
  MM1 = M1[0];
  MM2 = M2[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  ixmin2 = 0;
  ixmax2 = MM2 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
  xdelta2 = (b2[0] - a2[0]) / (MM2 - 1);
 
  // set all est = 0 
  for (i=0; i < MM1*MM2; i++)
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i]) && R_FINITE(x2[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      xpos2 = (x2[i] - a2[0]) / xdelta2;
      ix1 = floor(xpos1);
      ix2 = floor(xpos2);
      fx1 = xpos1 - ix1;
      fx2 = xpos2 - ix2;
      wi = weight[i];   
      
      if(ixmin1 <= ix1 && ixmin2 <= ix2 && ix2 <= ixmax2) {
	est[ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2);   
	est[ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2);
	est[(ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2;   
	est[(ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2;
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2) {
        est[ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2);   
	est[(ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2;   
      }
      else if (ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1) {
	est[ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2);   
	est[ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2);
      }
      else if (ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1) {
	est[ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2);  
      }
   }
  } 
}

void massdist3d(double *x1, double *x2,	double *x3, int *n, 
		double *a1, double *a2, double *a3, 
		double *b1, double *b2, double *b3,
		int *M1, int *M2, int *M3, double *weight, double *est)
{
  double fx1, fx2, fx3, xdelta1, xdelta2, xdelta3, xpos1, xpos2, xpos3, wi;   
  int i, ix1, ix2, ix3, ixmax1, ixmin1, ixmax2, ixmax3, ixmin2, ixmin3, MM1, MM2, MM3;
  
  MM1 = M1[0];
  MM2 = M2[0];
  MM3 = M3[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  ixmin2 = 0;
  ixmax2 = MM2 - 2;
  ixmin3 = 0;
  ixmax3 = MM3 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
  xdelta2 = (b2[0] - a2[0]) / (MM2 - 1);
  xdelta3 = (b3[0] - a3[0]) / (MM3 - 1);
 
  // set all est = 0 
  for (i=0; i < MM1*MM2*MM3; i++)  
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i]) && R_FINITE(x2[i]) && R_FINITE(x3[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      xpos2 = (x2[i] - a2[0]) / xdelta2;
      xpos3 = (x3[i] - a3[0]) / xdelta3;
      ix1 = floor(xpos1);
      ix2 = floor(xpos2);
      ix3 = floor(xpos3);
      fx1 = xpos1 - ix1;
      fx2 = xpos2 - ix2;
      fx3 = xpos3 - ix3;
      wi = weight[i];   
      
      if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
	est[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3);
	est[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3);   
	est[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3);
	est[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3;   
	est[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3;
	est[(ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3;   
	est[(ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*fx3;
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
	est[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3);   
	est[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3;   
	est[(ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3;   
      }   
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
	est[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3);
	est[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3;   
	est[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3;
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
	est[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3);
	est[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3);   
	est[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3);
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
 	est[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3;   
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
	est[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3);   
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   
	est[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3);
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1) {
	est[ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3);   	
      }
    }
  }
}

void massdist4d(double *x1, double *x2,	double *x3, double *x4, int *n, 
		double *a1, double *a2, double *a3, double *a4, 
		double *b1, double *b2, double *b3, double *b4, 
		int *M1, int *M2, int *M3, int *M4, double *weight, double *est)
{
  double fx1, fx2, fx3, fx4, xdelta1, xdelta2, xdelta3, xdelta4, xpos1, xpos2, xpos3, xpos4, wi;   
  int i, ix1, ix2, ix3, ix4, ixmax1, ixmin1, ixmax2, ixmax3, ixmax4, ixmin2, ixmin3, ixmin4, MM1, MM2, MM3, MM4;
  
  MM1 = M1[0];
  MM2 = M2[0];
  MM3 = M3[0];
  MM4 = M4[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  ixmin2 = 0;
  ixmax2 = MM2 - 2;
  ixmin3 = 0;
  ixmax3 = MM3 - 2;
  ixmin4 = 0;
  ixmax4 = MM4 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
  xdelta2 = (b2[0] - a2[0]) / (MM2 - 1);
  xdelta3 = (b3[0] - a3[0]) / (MM3 - 1);
  xdelta4 = (b4[0] - a4[0]) / (MM4 - 1);
 
  // set all est = 0 
  for (i=0; i < MM1*MM2*MM3*MM4; i++)  
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i]) && R_FINITE(x2[i]) && R_FINITE(x3[i]) && R_FINITE(x4[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      xpos2 = (x2[i] - a2[0]) / xdelta2;
      xpos3 = (x3[i] - a3[0]) / xdelta3;
      xpos4 = (x4[i] - a4[0]) / xdelta4;
      ix1 = floor(xpos1);
      ix2 = floor(xpos2);
      ix3 = floor(xpos3);
      ix4 = floor(xpos4);
      fx1 = xpos1 - ix1;
      fx2 = xpos2 - ix2;
      fx3 = xpos3 - ix3;
      fx4 = xpos4 - ix4;
      wi = weight[i];   
      
      if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3*(1-fx4);
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*fx3*(1-fx4);
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*fx4;
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3)*fx4;
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3*fx4;
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*fx3*fx4;
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4); 
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4);  
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3*(1-fx4);  
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4; 
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*fx4;  
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*fx4; 
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3*fx4;   
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3*(1-fx4);
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*fx4;
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3*fx4;
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3)*(1-fx4);
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*fx4;
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3)*fx4;
      } 
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3*(1-fx4);
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*fx3*(1-fx4);
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4); 
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4); 
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4; 
	est[(ix4+1)*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*fx4; 
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4); 
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4; 
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*fx4;   
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4); 
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4);  
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*fx3*(1-fx4); 
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4;   
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*fx4;
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*fx3*(1-fx4);
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1] += wi*fx1*fx2*(1-fx3)*(1-fx4);
      }  
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1 && ixmin4 <= ix4 && ix4 <= ixmax4) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4); 
	est[(ix4+1)*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*fx4; 
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4); 
	est[ix4*MM1*MM2*MM3 + (ix3+1)*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*fx3*(1-fx4); 
      } 
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + (ix2+1)*MM1 + ix1] += wi*(1-fx1)*fx2*(1-fx3)*(1-fx4); 
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4);   
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1 + 1] += wi*fx1*(1-fx2)*(1-fx3)*(1-fx4);
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1 && ix4 == ixmax4 + 1) {
	est[ix4*MM1*MM2*MM3 + ix3*MM1*MM2 + ix2*MM1 + ix1] += wi*(1-fx1)*(1-fx2)*(1-fx3)*(1-fx4); 
      }    
    }
  }
}

void interp1d(double *x1, int *n, double *a1, double *b1, int *M1,
              double *fun, double *est)
{ 
  double fx1, xdelta1, xpos1;   
  int i, ix1, ixmax1, ixmin1, MM1;

  MM1 = M1[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
   
  // set all est = 0 
  for (i=0; i < n[0]; i++)
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      ix1 = floor(xpos1);
      fx1 = xpos1 - ix1;   
      
      if(ixmin1 <= ix1 && ix1 <= ixmax1) {
	est[i] = fun[ix1]*(1-fx1) + fun[ix1 + 1]*fx1;
      }
      else if(ix1 <= -1) {
	est[i] = fun[0]; 
      }
      else if(ix1 >= ixmax1 + 1) {
	est[i] = fun[ixmax1 + 1];  
      }
    }
  } 
}

void interp2d(double *x1, double *x2,	int *n, 
	      double *a1, double *a2, double *b1, double *b2,
	      int *M1, int *M2, double *fun, double *est)
{
  double fx1, fx2, xdelta1, xdelta2, xpos1, xpos2;   
  int i, ix1, ix2, ixmax1, ixmin1, ixmax2, ixmin2, MM1, MM2;
  
  MM1 = M1[0];
  MM2 = M2[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  ixmin2 = 0;
  ixmax2 = MM2 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
  xdelta2 = (b2[0] - a2[0]) / (MM2 - 1);
 
  // set all est = 0 
  for (i=0; i < n[0]; i++)
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i]) && R_FINITE(x2[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      xpos2 = (x2[i] - a2[0]) / xdelta2;
      ix1 = floor(xpos1);
      ix2 = floor(xpos2);
      fx1 = xpos1 - ix1;
      fx2 = xpos2 - ix2;
           
      if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2) {
	est[i] = fun[ix2*MM1 + ix1]*(1-fx1)*(1-fx2) \
               + fun[ix2*MM1 + ix1 + 1]*fx1*(1-fx2) \
               + fun[(ix2+1)*MM1 + ix1]*(1-fx1)*fx2 \
               + fun[(ix2+1)*MM1 + ix1 + 1]*fx1*fx2;
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2) {
        est[i] = fun[ix2*MM1 + ix1]*(1-fx1)*(1-fx2) \
               + fun[(ix2+1)*MM1 + ix1]*(1-fx1)*fx2;  
      }
      else if (ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1) {
	est[i] = fun[ix2*MM1 + ix1]*(1-fx1)*(1-fx2) \
               + fun[ix2*MM1 + ix1 + 1]*fx1*(1-fx2);
      }
      else if (ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1) {
        est[i] = fun[ix2*MM1 + ix1]*(1-fx1)*(1-fx2);
      } 
    }
  } 
}

void interp3d(double *x1, double *x2,	double *x3, int *n, 
	      double *a1, double *a2, double *a3, 
	      double *b1, double *b2, double *b3,
              int *M1, int *M2, int *M3, double *fun, double *est)
{
  double fx1, fx2, fx3, xdelta1, xdelta2, xdelta3, xpos1, xpos2, xpos3;   
  int i, ix1, ix2, ix3, ixmax1, ixmin1, ixmax2, ixmax3, ixmin2, ixmin3, MM1, MM2, MM3;
  
  MM1 = M1[0];
  MM2 = M2[0];
  MM3 = M3[0];
  ixmin1 = 0;
  ixmax1 = MM1 - 2;
  ixmin2 = 0;
  ixmax2 = MM2 - 2;
  ixmin3 = 0;
  ixmax3 = MM3 - 2;
  xdelta1 = (b1[0] - a1[0]) / (MM1 - 1);
  xdelta2 = (b2[0] - a2[0]) / (MM2 - 1);
  xdelta3 = (b3[0] - a3[0]) / (MM3 - 1);
 
  // set all est = 0 
  for (i=0; i < n[0]; i++)  
    est[i] = 0.0;

  // assign linear binning weights
  for (i=0; i < n[0]; i++) {
    if(R_FINITE(x1[i]) && R_FINITE(x2[i]) && R_FINITE(x3[i])) {
      xpos1 = (x1[i] - a1[0]) / xdelta1;
      xpos2 = (x2[i] - a2[0]) / xdelta2;
      xpos3 = (x3[i] - a3[0]) / xdelta3;
      ix1 = floor(xpos1);
      ix2 = floor(xpos2);
      ix3 = floor(xpos3);
      fx1 = xpos1 - ix1;
      fx2 = xpos2 - ix2;
      fx3 = xpos3 - ix3;
      
      if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
               + fun[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1]*fx1*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1]*(1-fx1)*fx2*(1-fx3) \
	       + fun[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1]*fx1*fx2*(1-fx3) \
	       + fun[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*fx3 \
               + fun[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1]*fx1*(1-fx2)*fx3 \
	       + fun[(ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1]*(1-fx1)*fx2*fx3 \
	       + fun[(ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1]*fx1*fx2*fx3;
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1]*(1-fx1)*fx2*(1-fx3) \
	       + fun[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*fx3 \
	       + fun[(ix3+1)*MM1*MM2 + (ix2+1)*MM1 + ix1]*(1-fx1)*fx2*fx3;   
      }   
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1]*fx1*(1-fx2)*(1-fx3) \
	       + fun[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*fx3 \
	       + fun[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1 + 1]*fx1*(1-fx2)*fx3;
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1]*fx1*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1]*(1-fx1)*fx2*(1-fx3) \
	       + fun[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1 + 1]*fx1*fx2*(1-fx3);
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ixmin3 <= ix3 && ix3 <= ixmax3) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
 	       + fun[(ix3+1)*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*fx3;   
      }
      else if(ix1 == ixmax1 + 1 && ixmin2 <= ix2 && ix2 <= ixmax2 && ix3 == ixmax3 + 1) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + (ix2+1)*MM1 + ix1]*(1-fx1)*fx2*(1-fx3);   
      }
      else if(ixmin1 <= ix1 && ix1 <= ixmax1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3) \
	       + fun[ix3*MM1*MM2 + ix2*MM1 + ix1 + 1]*fx1*(1-fx2)*(1-fx3);
      }
      else if(ix1 == ixmax1 + 1 && ix2 == ixmax2 + 1 && ix3 == ixmax3 + 1) {
	est[i] = fun[ix3*MM1*MM2 + ix2*MM1 + ix1]*(1-fx1)*(1-fx2)*(1-fx3);   	
      }
    }
  }
}


/* Registration of native routines added 17/03/2017 */

static R_NativePrimitiveArgType md1_t[] = {
  REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType md2_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType md3_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType md4_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

const static R_CMethodDef cMethods[] = {
  {"massdist1d", (DL_FUNC) &massdist1d, 7,  md1_t},
  {"massdist2d", (DL_FUNC) &massdist2d, 11, md2_t},
  {"massdist3d", (DL_FUNC) &massdist3d, 15, md3_t},
  {"massdist4d", (DL_FUNC) &massdist4d, 19, md4_t},
  {"interp1d", (DL_FUNC) &interp1d, 7,  md1_t},
  {"interp2d", (DL_FUNC) &interp2d, 11, md2_t},
  {"interp3d", (DL_FUNC) &interp3d, 15, md3_t},
  {NULL, NULL, 0}
};


void attribute_visible R_init_ks(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE); 	
}
