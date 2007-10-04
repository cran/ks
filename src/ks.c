#include <stdlib.h>
#include <math.h>
#define E  2.7182818284590452354


double mult(double *x, double *vecA, double *y, int d);

/* 2-dim */
void dmvnorm_2d(double *x, double *y, double *mu, double *visigma, 
		double *detsigma, int *n, double *dens);
void dmvnormd1_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt);
void dmvnormd2_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt);
void dmvnormd3_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt);
void dmvnormd4_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt);
void dmvnormd5_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt);
void dmvnormd6_2d(double *x1, double *x2, double *vsigma, int *r, int *n,
		  double *derivt);

void dmvnorm_2d_sum(double *x1, double *x2, double *visigma,
		    double *detsigma, int *n, double *sum);
void dmvnorm_2d_sum_ab(double *x1, double *x2, double *sigma1, double *sigma2, 
					   int *n, int *which, double *sum);
void dmvnorm_2d_sum_clust(double *x1, double *x2, double *y1, double *y2, 
			  double *visigma, double *detsigma, int *nx, int *ny, 
			  double *sum);
void dmvnorm_2d_xxt_sum(double *x1, double *x2, double *visigma, 
			double *detsigma, int *n, double *sum);
void dmvnormd1_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);
void dmvnormd2_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);
void dmvnormd3_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);
void dmvnormd4_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);
void dmvnormd4_2d_xxt_sum(double *x1, double *x2, double *visigma, 
			  int*r, int *n, double *sum);
void dmvnormd5_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);
void dmvnormd6_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);

void dmvt_2d(double *x, double *y, double *mu, double *viSigma,  
	     double *df, int *n, double *dens);

/* 3-dim */
void dmvnorm_3d(double *x1, double *x2, double *x3, double *mu, 
		double *visigma, double *detsigma, int *n, double *dens);
void dmvnorm_3d_sum(double *x1, double *x2, double *x3, double *visigma, 
		    double *detsigma, int *n, double *sum);

/* 4-dim */
void dmvnorm_4d(double *x1, double *x2, double *x3, double *x4, double *mu, 
		double *visigma, double *detsigma, int *n, double *dens);
void dmvnormd2_4d(double *x1, double *x2, double *x3, double *x4, double *vsigma, 
		  int *r, int *n, double *derivt);
void dmvnorm_4d_sum(double *x1, double *x2, double *x3, double *x4,double *visigma, 
		    double *detsigma, int *n, double *sum);

/* 5-dim */
void dmvnorm_5d(double *x1, double *x2, double *x3, double *x4, double *x5, 
		double *mu, double *visigma, double *detsigma, 
		int *n, double *dens);
void dmvnorm_5d_sum(double *x1, double *x2, double *x3, double *x4, double *x5, 
		    double *visigma, double *detsigma, int *n, double *sum);

/* 6-dim */
void dmvnorm_6d(double *x1, double *x2, double *x3, double *x4, double *x5, 
		double *x6, double *mu, double *visigma, double *detsigma, 
		int *n, double *dens);
void dmvnorm_6d_sum(double *x1, double *x2, double *x3, double *x4, double *x5, 
		    double *x6, double *visigma, 
		    double *detsigma, int *n, double *sum);

/*************************************************************************
* Quadratic form
*
* Parameters
* x - x vector
* vecA - vector form of matrix A
* y - y vector
* d - dimension of X and Y
*
* Returns
* Quadratic form x^T * A * y
*************************************************************************/
double mult(double *x, double *vecA, double *y, int d)
{
    int i, j;
    double result;

    result = 0.0;
    for (j = 1; j <= d; j++)
      for (i = 1; i <= d; i++)
	result = result + x[i - 1] * vecA[(i - 1) + (j - 1)*d] * y[j - 1];
    
    return(result);
}


/*************************************************************************
 * Bivariate normal density
 *
 * Parameters
 * x - x values
 * y - y values
 * mu - mean
 * visigms - vector form of inverse of covariance matrix
 * detsigma - determinant of covariance matrix
 * n - number of values = length(X) = length(Y)
 * dens - contains density values
 *************************************************************************/
void dmvnorm_2d(double *x, double *y, double *mu, double *visigma, 
		double *detsigma, int *n, double *dens)
{
  double norm;
  double xmu[2];
  int j, d;
  
  d = 2;
  norm = 1/sqrt(pow(2*M_PI, d) * detsigma[0]);
  for (j = 1; j <= n[0]; j++)
  {
    xmu[0] = x[j - 1] - mu[0];
    xmu[1] = y[j - 1] - mu[1];
    dens[j - 1] = norm * exp(-0.5 * mult(xmu, visigma, xmu, d));
  }
}


/*************************************************************************
* Bivariate normal density - 1st order derivatives
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* x - number of values
* derivt - contains density derivative values
*************************************************************************/
void dmvnormd1_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt)
{
  double sigma11, sigma22, rho, x, y;
  int i;
  
  sigma11 = sqrt(vsigma[0]);
  sigma22 = sqrt(vsigma[3]);
  rho = vsigma[1] /(sigma11 * sigma22);


  /* first order derivatives */
  
  if ((r[0] == 1) && (r[1] == 0)) 
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      
      derivt[i - 1] = 
	(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
         pow(sigma11,2)*pow(y,2))/
       (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
     (2*pow(sigma22,2)*x - 2*rho*sigma11*sigma22*y))/
      (4.*M_PI*sqrt(1 - pow(rho,2))*(-1 + pow(rho,2))*pow(sigma11,3)*pow(sigma22,3));
    }
  else if ((r[0] == 0) && (r[1] == 1))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      
      derivt[i - 1] =
	(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
         pow(sigma11,2)*pow(y,2))/
       (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
     (-2*rho*sigma11*sigma22*x + 2*pow(sigma11,2)*y))/
      (4.*M_PI*sqrt(1 - pow(rho,2))*(-1 + pow(rho,2))*pow(sigma11,3)*pow(sigma22,3));
    }
}

/*************************************************************************
* Bivariate normal density - 2nd order derivatives
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* x - number of values
* derivt - contains density derivative values
*************************************************************************/
void dmvnormd2_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt)
{
    double sigma11, sigma22, rho, x, y;
    int i;
    
    sigma11 = sqrt(vsigma[0]);
    sigma22 = sqrt(vsigma[3]);
    rho = vsigma[1] /(sigma11 * sigma22);
    
    /* second order derivatives */
    
    if ((r[0] == 2) && (r[1] == 0)) 
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] = (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (pow(sigma22,2)*pow(x,2) - 
	      2*rho*sigma11*sigma22*x*y + pow(sigma11,2)*
	      ((-1 + pow(rho,2))*pow(sigma22,2) + 
	       pow(rho,2)*pow(y,2))))/
	    (2.*M_PI*pow(1 - pow(rho,2),2.5)*pow(sigma11,5)*
	     pow(sigma22,3));
	}
    else if ((r[0] == 1) && (r[1] == 1))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] = 
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (-(pow(rho,3)*pow(sigma11,2)*pow(sigma22,2)) + 
	      sigma11*sigma22*x*y + pow(rho,2)*sigma11*sigma22*x*y + 
	      rho*(-(pow(sigma22,2)*pow(x,2)) + 
		   pow(sigma11,2)*(pow(sigma22,2) - pow(y,2))))
	     )/(2.*M_PI*pow(1 - pow(rho,2),2.5)*pow(sigma11,4)*
		pow(sigma22,4));
	}
    else if ((r[0] == 0) && (r[1] == 2))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] = 
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (pow(rho,2)*pow(sigma22,2)*pow(x,2) - 
	      2*rho*sigma11*sigma22*x*y + pow(sigma11,2)*
	      ((-1 + pow(rho,2))*pow(sigma22,2) + pow(y,2)))
	     )/(2.*M_PI*pow(1 - pow(rho,2),2.5)*pow(sigma11,3)*
		pow(sigma22,5));
	}
    
}


/*************************************************************************
* Bivariate normal density - 3rd order derivatives
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* x - number of values
* derivt - contains density derivative values
*************************************************************************/
void dmvnormd3_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt)
{
  double sigma11, sigma22, rho, x, y;
  int i;
  
  sigma11 = sqrt(vsigma[0]);
  sigma22 = sqrt(vsigma[3]);
  rho = vsigma[1] /(sigma11 * sigma22);
  
  /* third order derivatives */
  
  if ((r[0] == 3) && (r[1] == 0))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
	
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (pow(sigma22,3)*pow(x,3) - 3*rho*sigma11*pow(sigma22,2)*pow(x,2)*y + 
        3*pow(sigma11,2)*sigma22*x*((-1 + pow(rho,2))*pow(sigma22,2) + 
           pow(rho,2)*pow(y,2)) - 
        rho*pow(sigma11,3)*y*(3*(-1 + pow(rho,2))*pow(sigma22,2) + 
           pow(rho,2)*pow(y,2))))/
	(2.*M_PI*pow(1 - pow(rho,2),3.5)*pow(sigma11,7)*pow(sigma22,4));
    }
  else if ((r[0] == 2) && (r[1] == 1))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
	
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (2*pow(rho,4)*pow(sigma11,3)*pow(sigma22,2)*y + 
        sigma11*pow(sigma22,2)*(-pow(sigma11,2) + pow(x,2))*y - 
        pow(rho,3)*pow(sigma11,2)*sigma22*x*(3*pow(sigma22,2) + pow(y,2)) + 
        rho*sigma22*x*(-(pow(sigma22,2)*pow(x,2)) + 
           pow(sigma11,2)*(3*pow(sigma22,2) - 2*pow(y,2))) + 
        pow(rho,2)*(2*sigma11*pow(sigma22,2)*pow(x,2)*y + 
           pow(sigma11,3)*(-(pow(sigma22,2)*y) + pow(y,3)))))/
	   (2.*M_PI*pow(1 - pow(rho,2),3.5)*pow(sigma11,6)*pow(sigma22,5));
    }
  else if ((r[0] == 1) && (r[1] == 2))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
	
      derivt[i - 1] = (pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (pow(rho,2)*pow(sigma22,3)*pow(x,3) - 
        rho*(2 + pow(rho,2))*sigma11*pow(sigma22,2)*pow(x,2)*y + 
        (1 + 2*pow(rho,2))*pow(sigma11,2)*sigma22*x*
         ((-1 + pow(rho,2))*pow(sigma22,2) + pow(y,2)) - 
        rho*pow(sigma11,3)*y*(3*(-1 + pow(rho,2))*pow(sigma22,2) + pow(y,2))))/
      (2.*M_PI*pow(1 - pow(rho,2),3.5)*pow(sigma11,5)*pow(sigma22,6));
    }
  else if ((r[0] == 0) && (r[1] == 3))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
	
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (-(pow(rho,3)*pow(sigma22,3)*x*(3*pow(sigma11,2) + pow(x,2))) + 
        3*pow(rho,2)*sigma11*pow(sigma22,2)*(pow(sigma11,2) + pow(x,2))*y + 
        3*rho*pow(sigma11,2)*sigma22*x*(pow(sigma22,2) - pow(y,2)) + 
        pow(sigma11,3)*y*(-3*pow(sigma22,2) + pow(y,2))))/
      (2.*M_PI*pow(1 - pow(rho,2),3.5)*pow(sigma11,4)*pow(sigma22,7));
    }
}


/*************************************************************************
* Bivariate normal density - 4th order derivatives
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* x - number of values
* derivt - contains density derivative values
*************************************************************************/
void dmvnormd4_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt)
{
  double sigma11, sigma22, rho, x, y;
  int i;
  
  sigma11 = sqrt(vsigma[0]);
  sigma22 = sqrt(vsigma[3]);
  rho = vsigma[1] /(sigma11 * sigma22);
  
  /* fourth order derivatives */
  
  if ((r[0] == 4) && (r[1] == 0))
    for (i = 1; i <= n[0]; i++)
      {
	x = x1[i - 1];
	y = x2[i - 1];
	
	derivt[i - 1] = 
	  (pow(M_E,(pow(sigma22,2)*pow(x,2) -2*rho*sigma11*sigma22*x*y + 
		    pow(sigma11,2)*pow(y,2))/
	       (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		pow(sigma22,2)))*
	   (pow(sigma22,4)*pow(x,4) - 
	    4*rho*sigma11*pow(sigma22,3)*pow(x,3)*y + 
	    6*pow(sigma11,2)*pow(sigma22,2)*pow(x,2)*
	    ((-1 + pow(rho,2))*pow(sigma22,2) + 
	     pow(rho,2)*pow(y,2)) - 4*rho*pow(sigma11,3)*sigma22*x*y*
	    (3*(-1 + pow(rho,2))*pow(sigma22,2) + 
	     pow(rho,2)*pow(y,2)) + pow(sigma11,4)*
	    (3*pow(-1 + pow(rho,2),2)*pow(sigma22,4) + 
	     6*pow(rho,2)*(-1 + pow(rho,2))*pow(sigma22,2)*pow(y,2) + 
	     pow(rho,4)*pow(y,4))))/
	  (2.*M_PI*pow(1 - pow(rho,2),4.5)*pow(sigma11,9)*
	       pow(sigma22,5));
      }

  else if ((r[0] == 3) && (r[1] == 1))
    for (i = 1; i <= n[0]; i++)
      {
	x = x1[i - 1];
	y = x2[i - 1];
	
	derivt[i - 1] =  
	  (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		    2*rho*sigma11*sigma22*x*y + 
		    pow(sigma11,2)*pow(y,2))/
		     (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		      pow(sigma22,2)))*
	   (sigma11*pow(sigma22,3)*x*
	    (-3*pow(sigma11,2) + pow(x,2))*y - 
	    3*pow(rho,5)*pow(sigma11,4)*pow(sigma22,2)*
	    (pow(sigma22,2) + pow(y,2)) + 
	    pow(rho,4)*pow(sigma11,3)*sigma22*x*y*
	    (9*pow(sigma22,2) + pow(y,2)) + 
	    3*pow(rho,2)*sigma11*sigma22*x*y*
	    (pow(sigma22,2)*pow(x,2) +pow(sigma11,2)*
	     (-2*pow(sigma22,2) + pow(y,2))) - 
	    rho*pow(sigma22,2)*
	    (pow(sigma22,2)*pow(x,4) + 3*pow(sigma11,4)*
	     (pow(sigma22,2) - pow(y,2)) + 
	     3*pow(sigma11,2)*pow(x,2)*
	     (-2*pow(sigma22,2) + pow(y,2))) + 
	    pow(rho,3)*pow(sigma11,2)*
	    (-3*pow(sigma22,2)*pow(x,2)*
	     (2*pow(sigma22,2) + pow(y,2)) + 
	     pow(sigma11,2)*(6*pow(sigma22,4) - pow(y,4))
	     )))/(2.*M_PI*pow(1 - pow(rho,2),4.5)*pow(sigma11,8)*
		  pow(sigma22,6));
      }
  
  else if ((r[0] == 2) && (r[1] == 2))
    for (i = 1; i <= n[0]; i++)
      {
            x = x1[i - 1];
            y = x2[i - 1];
	    derivt[i - 1] =
	      (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
			2*rho*sigma11*sigma22*x*y + 
			pow(sigma11,2)*pow(y,2))/
		   (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		    pow(sigma22,2)))*
	       (pow(rho,2)*pow(sigma22,4)*pow(x,4) - 
		2*rho*(1 + pow(rho,2))*sigma11*pow(sigma22,3)*
		pow(x,3)*y - 2*rho*pow(sigma11,3)*sigma22*x*y*
		(2*(-2 + pow(rho,2) + pow(rho,4))*
		 pow(sigma22,2) + (1 + pow(rho,2))*pow(y,2))
		+ pow(sigma11,2)*pow(sigma22,2)*pow(x,2)*
		((-1 - 4*pow(rho,2) + 5*pow(rho,4))*
		 pow(sigma22,2) + 
		 (1 + 4*pow(rho,2) + pow(rho,4))*pow(y,2)) + 
		pow(sigma11,4)*
		((1 - 3*pow(rho,4) + 2*pow(rho,6))*
		 pow(sigma22,4) + 
		 (-1 - 4*pow(rho,2) + 5*pow(rho,4))*
		 pow(sigma22,2)*pow(y,2) + 
	 pow(rho,2)*pow(y,4))))/
	      (2.*M_PI*pow(1 - pow(rho,2),4.5)*pow(sigma11,7)*
	       pow(sigma22,7));
      }

  else if ((r[0] == 1) && (r[1] == 3))
    for (i = 1; i <= n[0]; i++)
      {
	x = x1[i - 1];
	y = x2[i - 1];
	
	derivt[i - 1] = 
	  (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		    2*rho*sigma11*sigma22*x*y + 
		    pow(sigma11,2)*pow(y,2))/
	       (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		pow(sigma22,2)))*
	   (-3*pow(rho,5)*pow(sigma11,2)*pow(sigma22,4)*
	    (pow(sigma11,2) + pow(x,2)) + 
	    pow(rho,4)*sigma11*pow(sigma22,3)*x*
	    (9*pow(sigma11,2) + pow(x,2))*y + 
	    pow(sigma11,3)*sigma22*x*y*
	    (-3*pow(sigma22,2) + pow(y,2)) + 
	    pow(rho,3)*pow(sigma22,2)*
	    (-(pow(sigma22,2)*pow(x,4)) - 
	     3*pow(sigma11,2)*pow(x,2)*pow(y,2) + 
	     6*pow(sigma11,4)*(pow(sigma22,2) - pow(y,2))
	     ) + 3*pow(rho,2)*sigma11*sigma22*x*y*
	    (pow(sigma22,2)*pow(x,2) + 
	     pow(sigma11,2)*
	     (-2*pow(sigma22,2) + pow(y,2))) - 
	    rho*pow(sigma11,2)*
	    (3*pow(sigma22,2)*pow(x,2)*
	     (-pow(sigma22,2) + pow(y,2)) + 
	     pow(sigma11,2)*
	     (3*pow(sigma22,4) - 
	      6*pow(sigma22,2)*pow(y,2) + pow(y,4)))))/
	  (2.*M_PI*pow(1 - pow(rho,2),4.5)*pow(sigma11,6)*
	   pow(sigma22,8));
      }
  
    else if ((r[0] == 0) && (r[1] == 4))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] =
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		   (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		    pow(sigma22,2)))*
	     (pow(rho,4)*pow(sigma22,4)*pow(x,4) - 
	      4*pow(rho,3)*sigma11*pow(sigma22,3)*pow(x,3)*
	      y + 6*pow(rho,2)*pow(sigma11,2)*
	      pow(sigma22,2)*pow(x,2)*
	      ((-1 + pow(rho,2))*pow(sigma22,2) + pow(y,2))
	      - 4*rho*pow(sigma11,3)*sigma22*x*y*
	      (3*(-1 + pow(rho,2))*pow(sigma22,2) +pow(y,2)) + 
	      pow(sigma11,4)*
		(3*pow(-1 + pow(rho,2),2)*pow(sigma22,4) + 
		 6*(-1 + pow(rho,2))*pow(sigma22,2)*
		 pow(y,2) + pow(y,4))))/
	    (2.*M_PI*pow(1 - pow(rho,2),4.5)*pow(sigma11,5)*pow(sigma22,9));
	}
}


/*************************************************************************
* Bivariate normal density - 5th order derivatives
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* x - number of values
* derivt - contains density derivative values
*************************************************************************/
void dmvnormd5_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt)
{
  /* fifth order derivatives */

  double sigma11, sigma22, rho, x, y;
  int i;
  
  sigma11 = sqrt(vsigma[0]);
  sigma22 = sqrt(vsigma[3]);
  rho = vsigma[1] /(sigma11 * sigma22);
    
  if ((r[0] == 5) && (r[1] == 0))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (sigma22*x - rho*sigma11*y)*(240*pow(-1 + pow(rho,2),2)*pow(sigma11,4)*
         pow(sigma22,8) + 160*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,6)*
         pow(sigma22*x - rho*sigma11*y,2) + 
        pow(2*pow(sigma22,2)*x - 2*rho*sigma11*sigma22*y,4)))/
	(32.*M_PI*pow(1 - pow(rho,2),5.5)*pow(sigma11,11)*pow(sigma22,10));
    }
  else if ((r[0] == 4) && (r[1] == 1))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (-3*pow(-1 + pow(rho,2),2)*pow(sigma11,4)*pow(sigma22,4)*
         (rho*sigma22*x - sigma11*y) - 
        6*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)*(rho*sigma22*x - sigma11*y)*
         pow(sigma22*x - rho*sigma11*y,2) + 
        (-(rho*sigma22*x) + sigma11*y)*pow(sigma22*x - rho*sigma11*y,4) + 
        12*rho*pow(-1 + pow(rho,2),2)*pow(sigma11,4)*pow(sigma22,4)*
         (-(sigma22*x) + rho*sigma11*y) + 
        4*rho*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)*
         pow(-(sigma22*x) + rho*sigma11*y,3)))/
	(2.*M_PI*pow(1 - pow(rho,2),5.5)*pow(sigma11,10)*pow(sigma22,7));
    }
  else if ((r[0] == 3) && (r[1] == 2))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (pow(rho,2)*pow(sigma22,5)*pow(x,5) - 
        rho*(2 + 3*pow(rho,2))*sigma11*pow(sigma22,4)*pow(x,4)*y - 
        rho*pow(sigma11,3)*pow(sigma22,2)*pow(x,2)*y*
         (15*(-1 + pow(rho,4))*pow(sigma22,2) + 
           (3 + 6*pow(rho,2) + pow(rho,4))*pow(y,2)) + 
        pow(sigma11,2)*pow(sigma22,3)*pow(x,3)*
         ((-1 - 8*pow(rho,2) + 9*pow(rho,4))*pow(sigma22,2) + 
           (1 + 6*pow(rho,2) + 3*pow(rho,4))*pow(y,2)) - 
        rho*pow(sigma11,5)*y*(3*pow(-1 + pow(rho,2),2)*(3 + 2*pow(rho,2))*
            pow(sigma22,4) + (-3 - 4*pow(rho,2) + 7*pow(rho,4))*pow(sigma22,2)*
            pow(y,2) + pow(rho,2)*pow(y,4)) + 
        pow(sigma11,4)*sigma22*x*(3*pow(-1 + pow(rho,2),2)*(1 + 4*pow(rho,2))*
            pow(sigma22,4) + 3*(-1 - 6*pow(rho,2) + 5*pow(rho,4) + 2*pow(rho,6))*
            pow(sigma22,2)*pow(y,2) + pow(rho,2)*(3 + 2*pow(rho,2))*pow(y,4))))/
	(2.*M_PI*pow(1 - pow(rho,2),5.5)*pow(sigma11,9)*pow(sigma22,8));
    }
  else if ((r[0] == 2) && (r[1] == 3))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      derivt[i - 1] = (pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
         pow(sigma11,2)*pow(y,2))/
       (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
     (6*pow(rho,7)*pow(sigma11,4)*pow(sigma22,5)*x - 
       6*pow(rho,6)*pow(sigma11,3)*pow(sigma22,4)*(2*pow(sigma11,2) + pow(x,2))*
        y - pow(sigma11,3)*pow(sigma22,2)*(pow(sigma11,2) - pow(x,2))*y*
        (3*pow(sigma22,2) - pow(y,2)) + 
       pow(rho,4)*sigma11*pow(sigma22,2)*y*
        (-2*pow(sigma22,2)*pow(x,4) + 
          3*pow(sigma11,4)*(7*pow(sigma22,2) - 3*pow(y,2)) - 
          3*pow(sigma11,2)*pow(x,2)*(5*pow(sigma22,2) + pow(y,2))) + 
       pow(rho,5)*pow(sigma11,2)*pow(sigma22,3)*x*
        (-3*pow(sigma11,2)*(pow(sigma22,2) - 5*pow(y,2)) + 
          pow(x,2)*(7*pow(sigma22,2) + pow(y,2))) + 
       pow(rho,3)*(pow(sigma22,5)*pow(x,5) + 
          2*pow(sigma11,2)*pow(sigma22,3)*pow(x,3)*
           (-2*pow(sigma22,2) + 3*pow(y,2)) + 
          3*pow(sigma11,4)*sigma22*x*(-4*pow(sigma22,4) + pow(y,4))) - 
       pow(rho,2)*sigma11*y*(3*pow(sigma22,4)*pow(x,4) + 
          6*pow(sigma11,2)*pow(sigma22,2)*pow(x,2)*
           (-3*pow(sigma22,2) + pow(y,2)) + 
          pow(sigma11,4)*(6*pow(sigma22,4) - 8*pow(sigma22,2)*pow(y,2) + pow(y,4)))
         + rho*pow(sigma11,2)*sigma22*x*
        (3*pow(sigma22,2)*pow(x,2)*(-pow(sigma22,2) + pow(y,2)) + 
          pow(sigma11,2)*(9*pow(sigma22,4) - 15*pow(sigma22,2)*pow(y,2) + 
             2*pow(y,4)))))/
	(2.*M_PI*pow(1 - pow(rho,2),5.5)*pow(sigma11,8)*pow(sigma22,9));
    }
  else if ((r[0] == 1) && (r[1] == 4))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];  
      derivt[i - 1] = (pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
         pow(sigma11,2)*pow(y,2))/
       (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
     (-(pow(sigma11,2)*(sigma22*x - rho*sigma11*y)*
          (3*pow(-1 + pow(rho,2),2)*pow(sigma22,4) + 
            (6*(-1 + pow(rho,2))*pow(sigma22,2)*pow(rho*sigma22*x - sigma11*y,2))/
             pow(sigma11,2) + pow(rho*sigma22*x - sigma11*y,4)/pow(sigma11,4))) + 
       4*rho*(1 - pow(rho,2))*pow(sigma22,2)*
        (pow(rho,3)*pow(sigma22,3)*x*(3*pow(sigma11,2) + pow(x,2)) - 
          3*pow(rho,2)*sigma11*pow(sigma22,2)*(pow(sigma11,2) + pow(x,2))*y - 
          pow(sigma11,3)*y*(-3*pow(sigma22,2) + pow(y,2)) + 
          3*rho*pow(sigma11,2)*sigma22*x*(-pow(sigma22,2) + pow(y,2)))))/
	(2.*M_PI*pow(1 - pow(rho,2),5.5)*pow(sigma11,5)*pow(sigma22,10));
    }
  else if ((r[0] == 0) && (r[1] == 5))
    for (i = 1; i <= n[0]; i++)
    {
      x = x1[i - 1];
      y = x2[i - 1];
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
          pow(sigma11,2)*pow(y,2))/
        (2.*(-1 + pow(rho,2))*pow(sigma11,2)*pow(sigma22,2)))*
      (-(rho*sigma22*x) + sigma11*y)*(240*pow(-1 + pow(rho,2),2)*pow(sigma11,8)*
         pow(sigma22,4) + 160*(-1 + pow(rho,2))*pow(sigma11,6)*pow(sigma22,2)*
         pow(rho*sigma22*x - sigma11*y,2) + 
        pow(2*rho*sigma11*sigma22*x - 2*pow(sigma11,2)*y,4)))/
	(32.*M_PI*pow(1 - pow(rho,2),5.5)*pow(sigma11,10)*pow(sigma22,11));
    }
}

/*************************************************************************
* Bivariate normal density - 6th order derivatives
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* x - number of values
* derivt - contains density derivative values
*************************************************************************/
void dmvnormd6_2d(double *x1, double *x2, double *vsigma, int *r, int *n, 
		  double *derivt)
{
  /* sixth order derivatives */

  double sigma11, sigma22, rho, x, y;
  int i;

    sigma11 = sqrt(vsigma[0]);
    sigma22 = sqrt(vsigma[3]);
    rho = vsigma[1] /(sigma11 * sigma22);
    if ((r[0] == 6) && (r[1] == 0))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] = 
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) -2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)* pow(sigma22,2)))*
	     (pow(sigma22,6)*pow(x,6) - 
	      6*rho*sigma11*pow(sigma22,5)*pow(x,5)*y + 
	      15*pow(sigma11,2)*pow(sigma22,4)*pow(x,4)*
	      ((-1 + pow(rho,2))*pow(sigma22,2) + 
	       pow(rho,2)*pow(y,2)) - 
	      20*rho*pow(sigma11,3)*pow(sigma22,3)*pow(x,3)*
	      y*(3*(-1 + pow(rho,2))*pow(sigma22,2) + 
		 pow(rho,2)*pow(y,2)) + 
	      15*pow(sigma11,4)*pow(sigma22,2)*pow(x,2)*
	      (3*pow(-1 + pow(rho,2),2)*pow(sigma22,4) + 
	       6*pow(rho,2)*(-1 + pow(rho,2))*
	       pow(sigma22,2)*pow(y,2) + 
	       pow(rho,4)*pow(y,4)) - 
	      6*rho*pow(sigma11,5)*sigma22*x*y*
	      (15*pow(-1 + pow(rho,2),2)*pow(sigma22,4) + 
	       10*pow(rho,2)*(-1 + pow(rho,2))*
	       pow(sigma22,2)*pow(y,2) + 
	       pow(rho,4)*pow(y,4))+ pow(sigma11,6)*
	      (15*pow(-1 + pow(rho,2),3)*pow(sigma22,6) + 
	       45*pow(rho,2)*pow(-1 + pow(rho,2),2)*
	       pow(sigma22,4)*pow(y,2) + 
	       15*pow(rho,4)*(-1 + pow(rho,2))*
	       pow(sigma22,2)*pow(y,4) + 
	       pow(rho,6)*pow(y,6))))/
	    (2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,13)*
	     pow(sigma22,7));
        }
    
    else if ((r[0] == 5) && (r[1] == 1))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] =
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (sigma11*pow(sigma22,5)*x*
	      (15*pow(sigma11,4) - 
	       10*pow(sigma11,2)*pow(x,2) + pow(x,4))*y - 
	      5*pow(rho,7)*pow(sigma11,6)*pow(sigma22,2)*
	      (3*pow(sigma22,4) + 
	       6*pow(sigma22,2)*pow(y,2) + pow(y,4)) + 
	      pow(rho,6)*pow(sigma11,5)*sigma22*x*y*
	      (75*pow(sigma22,4) + 
	       30*pow(sigma22,2)*pow(y,2) + pow(y,4)) + 
	      5*pow(rho,2)*sigma11*pow(sigma22,3)*x*y*
	      (pow(sigma22,2)*pow(x,4) + 
	       pow(sigma11,4)*
	       (9*pow(sigma22,2) - 6*pow(y,2)) + 
	       2*pow(sigma11,2)*pow(x,2)*
	       (-4*pow(sigma22,2) + pow(y,2))) + 
	      rho*pow(sigma22,4)*
	      (-(pow(sigma22,2)*pow(x,6)) + 
	       15*pow(sigma11,6)*
	       (pow(sigma22,2) - pow(y,2)) + 
	       5*pow(sigma11,2)*pow(x,4)*
	       (3*pow(sigma22,2) - pow(y,2)) + 
	       15*pow(sigma11,4)*pow(x,2)*
	       (-3*pow(sigma22,2) + 2*pow(y,2))) - 
	      5*pow(rho,3)*pow(sigma11,2)*pow(sigma22,2)*
	      (pow(sigma22,2)*pow(x,4)*
	       (3*pow(sigma22,2) + 2*pow(y,2)) + 
	       pow(sigma11,4)*
	       (9*pow(sigma22,4) - 2*pow(y,4)) - 
	       2*pow(sigma11,2)*pow(x,2)*
	       (9*pow(sigma22,4) + 
		3*pow(sigma22,2)*pow(y,2) - pow(y,4))) + 
	      5*pow(rho,4)*pow(sigma11,3)*sigma22*x*y*
	      (2*pow(sigma22,2)*pow(x,2)*
	       (5*pow(sigma22,2) + pow(y,2)) + 
	       pow(sigma11,2)*
	       (-27*pow(sigma22,4) + pow(y,4))) + 
	      pow(rho,5)*pow(sigma11,4)*
	      (-5*pow(sigma22,2)*pow(x,2)*
	       (9*pow(sigma22,4) + 
		12*pow(sigma22,2)*pow(y,2) + pow(y,4)) + 
	       pow(sigma11,2)*
	       (45*pow(sigma22,6) + 
		45*pow(sigma22,4)*pow(y,2) - 
		5*pow(sigma22,2)*pow(y,4) - pow(y,6)))))/
	    (2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,12)*
	     pow(sigma22,8));
        }
    
    else if ((r[0] == 4) && (r[1] == 2))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  derivt[i - 1] =
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (pow(rho,2)*pow(sigma22,6)*pow(x,6) - 
	      2*rho*(1 + 2*pow(rho,2))*sigma11*pow(sigma22,5)*
	      pow(x,5)*y - 4*rho*pow(sigma11,3)*
	      pow(sigma22,3)*pow(x,3)*y*
	      ((-6 - 3*pow(rho,2) + 9*pow(rho,4))*
	       pow(sigma22,2) + 
	       (1 + 3*pow(rho,2) + pow(rho,4))*pow(y,2)) + 
	      pow(sigma11,2)*pow(sigma22,4)*pow(x,4)*
	      ((-1 - 13*pow(rho,2) + 14*pow(rho,4))*
	       pow(sigma22,2) + 
	       (1 + 8*pow(rho,2) + 6*pow(rho,4))*pow(y,2))
	      - 2*rho*pow(sigma11,5)*sigma22*x*y*
	      (3*pow(-1 + pow(rho,2),2)*(7 + 8*pow(rho,2))*
	       pow(sigma22,4) + 
	       2*(-3 - 7*pow(rho,2) + 8*pow(rho,4) + 
		  2*pow(rho,6))*pow(sigma22,2)*pow(y,2) + 
	       pow(rho,2)*(2 + pow(rho,2))*pow(y,4)) + 
	      pow(sigma11,4)*pow(sigma22,2)*pow(x,2)*
	      (3*pow(-1 + pow(rho,2),2)*(2 + 13*pow(rho,2))*
	       pow(sigma22,4) + 
	       6*(-1 - 8*pow(rho,2) + 4*pow(rho,4) + 
		  5*pow(rho,6))*pow(sigma22,2)*pow(y,2) + 
	       pow(rho,2)*(6 + 8*pow(rho,2) + pow(rho,4))*
	       pow(y,4)) + 
	      pow(sigma11,6)*
	      (3*pow(-1 + pow(rho,2),3)*(1 + 4*pow(rho,2))*
	       pow(sigma22,6) + 
	       3*pow(-1 + pow(rho,2),2)*
	       (1 + 10*pow(rho,2) + 4*pow(rho,4))*
	       pow(sigma22,4)*pow(y,2) + 
	       3*pow(rho,2)*
	       (-2 - pow(rho,2) + 3*pow(rho,4))*
	       pow(sigma22,2)*pow(y,4) + 
	       pow(rho,4)*pow(y,6))))/
	    (2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,11)*
	     pow(sigma22,9));
	}
    
    else if ((r[0] == 3) && (r[1] == 3))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  derivt[i - 1] = 
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (-6*pow(rho,9)*pow(sigma11,6)*pow(sigma22,6) + 
	      18*pow(rho,8)*pow(sigma11,5)*pow(sigma22,5)*x*
	      y + pow(sigma11,3)*pow(sigma22,3)*x*
	      (3*pow(sigma11,2) - pow(x,2))*y*
	      (3*pow(sigma22,2) - pow(y,2)) + 
	      9*pow(rho,7)*pow(sigma11,4)*pow(sigma22,4)*
	      (pow(sigma11,2)*
	       (pow(sigma22,2) - 3*pow(y,2)) - 
	       pow(x,2)*(3*pow(sigma22,2) + pow(y,2))) + 
	      pow(rho,6)*pow(sigma11,3)*pow(sigma22,3)*x*y*
	      (pow(x,2)*(21*pow(sigma22,2) + pow(y,2)) + 
	       3*pow(sigma11,2)*
	       (9*pow(sigma22,2) + 7*pow(y,2))) + 
	      3*pow(rho,5)*pow(sigma11,2)*pow(sigma22,2)*
        (-(pow(sigma22,2)*pow(x,4)*
	   (4*pow(sigma22,2) + pow(y,2))) + 
	 pow(sigma11,4)*
	 (3*pow(sigma22,4) + 
	  12*pow(sigma22,2)*pow(y,2) - 4*pow(y,4))
	 + pow(sigma11,2)*pow(x,2)*
	 (12*pow(sigma22,4) - 
	  15*pow(sigma22,2)*pow(y,2) - pow(y,4)))
	      + 3*pow(rho,2)*sigma11*sigma22*x*y*
	      (pow(sigma22,4)*pow(x,4) + 
	       pow(sigma11,2)*pow(sigma22,2)*pow(x,2)*
	       (-11*pow(sigma22,2) + 3*pow(y,2)) + 
	       pow(sigma11,4)*
	       (15*pow(sigma22,4) - 
		11*pow(sigma22,2)*pow(y,2) + pow(y,4)))
	      + 3*rho*pow(sigma11,2)*pow(sigma22,2)*
	      (pow(sigma22,2)*pow(x,4)*
	       (pow(sigma22,2) - pow(y,2)) - 
	       pow(sigma11,2)*pow(x,2)*
	       (6*pow(sigma22,4) - 
		9*pow(sigma22,2)*pow(y,2) + pow(y,4)) + 
	       pow(sigma11,4)*
	       (3*pow(sigma22,4) - 
		6*pow(sigma22,2)*pow(y,2) + pow(y,4))) + 
	      3*pow(rho,4)*sigma11*sigma22*x*y*
	      (pow(sigma22,4)*pow(x,4) + 
	       pow(sigma11,2)*pow(sigma22,2)*pow(x,2)*
	       (5*pow(sigma22,2) + 3*pow(y,2)) + 
	       pow(sigma11,4)*
	       (-33*pow(sigma22,4) + 
		5*pow(sigma22,2)*pow(y,2) + pow(y,4))) - 
	      pow(rho,3)*(pow(sigma22,6)*pow(x,6) + 
			  9*pow(sigma11,2)*pow(sigma22,4)*pow(x,4)*
			  (-pow(sigma22,2) + pow(y,2)) - 
			  9*pow(sigma11,4)*pow(sigma22,2)*pow(x,2)*
			  (pow(sigma22,4) + 
			   3*pow(sigma22,2)*pow(y,2) - pow(y,4)) + 
			  pow(sigma11,6)*
			  (21*pow(sigma22,6) - 
			   9*pow(sigma22,4)*pow(y,2) - 
			   9*pow(sigma22,2)*pow(y,4) + pow(y,6)))))/
	    (2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,10)*
	     pow(sigma22,10));
        }
    
    else if ((r[0] == 2) && (r[1] == 4))
      for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];  
	  derivt[i - 1] =
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (pow(rho,4)*pow(sigma22,6)*pow(x,6) - 
	      2*pow(rho,3)*(2 + pow(rho,2))*sigma11*
	      pow(sigma22,5)*pow(x,5)*y - 
	      4*rho*pow(sigma11,3)*pow(sigma22,3)*pow(x,3)*y*
	      ((-3 - 7*pow(rho,2) + 8*pow(rho,4) + 
		2*pow(rho,6))*pow(sigma22,2) + 
	       (1 + 3*pow(rho,2) + pow(rho,4))*pow(y,2)) + 
	      pow(rho,2)*pow(sigma11,2)*pow(sigma22,4)*
	      pow(x,4)*((-6 - 3*pow(rho,2) + 9*pow(rho,4))*
			pow(sigma22,2) + 
			(6 + 8*pow(rho,2) + pow(rho,4))*pow(y,2)) - 
	      2*rho*pow(sigma11,5)*sigma22*x*y*
	      (3*pow(-1 + pow(rho,2),2)*(7 + 8*pow(rho,2))*
	       pow(sigma22,4) + 
	       6*(-2 - pow(rho,2) + 3*pow(rho,4))*
	       pow(sigma22,2)*pow(y,2) + 
	       (1 + 2*pow(rho,2))*pow(y,4)) + 
	      pow(sigma11,4)*pow(sigma22,2)*pow(x,2)*
	      (3*pow(-1 + pow(rho,2),2)*
	       (1 + 10*pow(rho,2) + 4*pow(rho,4))*
	       pow(sigma22,4) + 
	       6*(-1 - 8*pow(rho,2) + 4*pow(rho,4) + 
		  5*pow(rho,6))*pow(sigma22,2)*pow(y,2) + 
	       (1 + 8*pow(rho,2) + 6*pow(rho,4))*pow(y,4))
	      + pow(sigma11,6)*
	      (3*pow(-1 + pow(rho,2),3)*(1 + 4*pow(rho,2))*
	       pow(sigma22,6) + 
	       3*pow(-1 + pow(rho,2),2)*
	       (2 + 13*pow(rho,2))*pow(sigma22,4)*pow(y,2)
	       + (-1 - 13*pow(rho,2) + 14*pow(rho,4))*
	       pow(sigma22,2)*pow(y,4) + 
	       pow(rho,2)*pow(y,6))))/
	    (2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,9)*
	     pow(sigma22,11));
        }

    else if ((r[0] == 1) && (r[1] == 5))
        for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  derivt[i - 1] =
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
		 (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
		  pow(sigma22,2)))*
	     (-5*pow(rho,7)*pow(sigma11,2)*pow(sigma22,6)*
        (3*pow(sigma11,4) + 
	 6*pow(sigma11,2)*pow(x,2) + pow(x,4)) + 
	      pow(rho,6)*sigma11*pow(sigma22,5)*x*
	      (75*pow(sigma11,4) + 
	       30*pow(sigma11,2)*pow(x,2) + pow(x,4))*y + 
	      pow(sigma11,5)*sigma22*x*y*
	      (15*pow(sigma22,4) - 
	       10*pow(sigma22,2)*pow(y,2) + pow(y,4)) + 
	      pow(rho,5)*pow(sigma22,4)*
	      (-(pow(sigma22,2)*pow(x,6)) + 
	       15*pow(sigma11,4)*pow(x,2)*
	       (3*pow(sigma22,2) - 4*pow(y,2)) + 
	       45*pow(sigma11,6)*
	       (pow(sigma22,2) - pow(y,2)) - 
	       5*pow(sigma11,2)*pow(x,4)*
	       (pow(sigma22,2) + pow(y,2))) + 
	      5*pow(rho,4)*sigma11*pow(sigma22,3)*x*y*
	      (pow(sigma22,2)*pow(x,4) + 
	       2*pow(sigma11,2)*pow(x,2)*pow(y,2) + 
	       pow(sigma11,4)*
	       (-27*pow(sigma22,2) + 10*pow(y,2))) + 
	      5*pow(rho,2)*pow(sigma11,3)*sigma22*x*y*
        (2*pow(sigma22,2)*pow(x,2)*
	 (-3*pow(sigma22,2) + pow(y,2)) + 
	 pow(sigma11,2)*
	 (9*pow(sigma22,4) - 
	  8*pow(sigma22,2)*pow(y,2) + pow(y,4))) - 
	      5*pow(rho,3)*pow(sigma11,2)*pow(sigma22,2)*
	      (2*pow(sigma11,2)*pow(x,2)*pow(y,2)*
	       (-3*pow(sigma22,2) + pow(y,2)) + 
	       2*pow(sigma22,2)*pow(x,4)*
	       (-pow(sigma22,2) + pow(y,2)) + 
	       3*pow(sigma11,4)*
	       (3*pow(sigma22,4) - 
		6*pow(sigma22,2)*pow(y,2) + pow(y,4))) + 
	      rho*pow(sigma11,4)*
	      (-5*pow(sigma22,2)*pow(x,2)*
	       (3*pow(sigma22,4) - 
		6*pow(sigma22,2)*pow(y,2) + pow(y,4)) + 
	       pow(sigma11,2)*
	       (15*pow(sigma22,6) - 
		45*pow(sigma22,4)*pow(y,2) + 
		15*pow(sigma22,2)*pow(y,4) - pow(y,6)))))
	    /(2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,8)*
	      pow(sigma22,12));
        }

    else if ((r[0] == 0) && (r[1] == 6))
        for (i = 1; i <= n[0]; i++)
        {
	  x = x1[i - 1];
	  y = x2[i - 1];
	  
	  derivt[i - 1] =
	    (pow(M_E,(pow(sigma22,2)*pow(x,2) - 
		      2*rho*sigma11*sigma22*x*y + 
		      pow(sigma11,2)*pow(y,2))/
     (2.*(-1 + pow(rho,2))*pow(sigma11,2)*
      pow(sigma22,2)))*
	     (pow(rho,6)*pow(sigma22,6)*pow(x,6) - 
	      6*pow(rho,5)*sigma11*pow(sigma22,5)*pow(x,5)*
	      y + 15*pow(rho,4)*pow(sigma11,2)*
	      pow(sigma22,4)*pow(x,4)*
	      ((-1 + pow(rho,2))*pow(sigma22,2) + pow(y,2))
	      - 20*pow(rho,3)*pow(sigma11,3)*pow(sigma22,3)*
	      pow(x,3)*y*(3*(-1 + pow(rho,2))*
			  pow(sigma22,2) + pow(y,2)) + 
	      15*pow(rho,2)*pow(sigma11,4)*pow(sigma22,2)*
	      pow(x,2)*(3*pow(-1 + pow(rho,2),2)*
			pow(sigma22,4) + 
			6*(-1 + pow(rho,2))*pow(sigma22,2)*
			pow(y,2) + pow(y,4)) - 
	      6*rho*pow(sigma11,5)*sigma22*x*y*
	      (15*pow(-1 + pow(rho,2),2)*pow(sigma22,4) + 
	       10*(-1 + pow(rho,2))*pow(sigma22,2)*
	       pow(y,2) + pow(y,4)) +pow(sigma11,6)*
	      (15*pow(-1 + pow(rho,2),3)*pow(sigma22,6) + 
	       45*pow(-1 + pow(rho,2),2)*pow(sigma22,4)*pow(y,2) + 
	       15*(-1 + pow(rho,2))*pow(sigma22,2)*
	       pow(y,4) + pow(y,6))))/
	    (2.*M_PI*pow(1 - pow(rho,2),6.5)*pow(sigma11,7)*
	     pow(sigma22,13));
        }
}

/*************************************************************************
* Double sum of normal density values - bivariate
*
* Parameters
* x1 - x values
* x2 - y values
* visigma - vec of inverse Sigma
* n - number of values
* sum - contains sum of density values
*************************************************************************/

void dmvnorm_2d_sum(double *x1, double *x2, double *visigma, double *detsigma,
    int *n, double *sum)
{
    int i, j, k, N;
    double mu[] = {0.0, 0.0};
    double sum1;
    double *x, *y, *dens;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    dens = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnorm_2d(x, y, mu, visigma, detsigma, n, dens);
        for (k = i; k <= N; k++)
            sum1 = sum1 + dens[k - 1];
    }

    sum[0] = sum1;
    free(x);
    free(y);
    free(dens);
}





/*************************************************************************
* Double sum of normal density values for Abramson's selector - bivariate
*
* Parameters
* x1 - x values
* x2 - y values
* visigma - vec of inverse Sigma
* n - number of values
* sum - contains sum of density values
*************************************************************************/

void dmvnorm_2d_sum_ab(double *x1, double *x2, double *sigma1, double *sigma2,
		       int *n, int *which, double *sum)
{
  int i, j, N, d, w;
  double visigma[] = {0.0, 0.0, 0.0, 0.0};
  double sum1, norm, detsigma;
  double *x;
  
  d = 2;
  N = n[0];
  w = which[0];
  sum1 = 0.0;
  x = malloc(sizeof(double)*d);

  for (i = 1; i <= N; i++)
  {
    for (j = 1; j <= N; j++)
    {
      x[0] = x1[i - 1] - x1[j - 1];
      x[1] = x2[i - 1] - x2[j - 1];
      if (w==1)
      {
	visigma[0] = 1/(sigma1[i - 1] + sigma1[j - 1]);
	visigma[3] = 1/(sigma2[i - 1] + sigma2[j - 1]);
      }
      else 
      {
	visigma[0] = 1/sigma1[j - 1];
	visigma[3] = 1/sigma2[j - 1];
      }
      detsigma = 1/(visigma[0]*visigma[3]); 
      norm = 1/sqrt(pow(2*M_PI, d) * detsigma); 
      if ((w==0) & (i==j))
	;
      else 
	sum1 = sum1 + norm * exp(-0.5 * mult(x, visigma, x, d));
    }
  }

  sum[0] = sum1;
  free(x);
}


/*************************************************************************
* Double sum of normal density values for re-clustered selector - bivariate
*
* Parameters
* x1 - x values
* x2 - y values
* visigma - vec of inverse Sigma
* n - number of values
* sum - contains sum of density values
*************************************************************************/

void dmvnorm_2d_sum_clust(double *x1, double *x2, double *y1, double *y2, 
			 double *visigma, double *detsigma, int *nx, int *ny, 
			 double *sum)
{
  int i, j, k, Nx, Ny;
  double mu[] = {0.0, 0.0};
  double sum1;
  double *x, *y, *dens;
    
  Nx = nx[0];
  Ny = ny[0];  
  x = malloc(sizeof(double)*Ny);
  y = malloc(sizeof(double)*Ny);
  dens = malloc(sizeof(double)*Ny);
  sum1 = 0.0;
  
  for (i = 1; i <= Nx; i++)
  {
    for (j = 1; j <= Ny; j++)
    {
      x[j - 1] = x1[i - 1] - y1[j - 1];
      y[j - 1] = x2[i - 1] - y2[j - 1];
    }
    dmvnorm_2d(x, y, mu, visigma, detsigma, ny, dens);
    for (k = 1; k <= Ny; k++)
      sum1 = sum1 + dens[k - 1];
  }

  sum[0] = sum1;
  free(x);
  free(y);
  free(dens);

}

/*************************************************************************
* Double sum of normal density values * x * x^T - bivariate
*
* Parameters
* x1 - x values
* x2 - y values
* visigma - vec of inverse Sigma
* detsigma - det of Sigma
* n - number of values
* sum - contains sum of density values
*************************************************************************/

void dmvnorm_2d_xxt_sum(double *x1, double *x2, double *visigma, 
     double *detsigma, int *n, double *sum)
{
    int i, j, k, N;
    double mu[] = {0.0, 0.0};
    double *Xmat0, *Xmat1, *Xmat2, *Xmat3;
    double *x, *y, *dens;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    dens = malloc(sizeof(double)*N);
    Xmat0 = malloc(sizeof(double)*N);
    Xmat1 = malloc(sizeof(double)*N);
    Xmat2 = malloc(sizeof(double)*N);
    Xmat3 = malloc(sizeof(double)*N);
     
    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
            Xmat0[j - 1] = x[j - 1] * x[j - 1];
            Xmat2[j - 1] = Xmat1[j - 1] = x[j - 1] * y[j - 1];
            Xmat3[j - 1] = y[j - 1] * y[j - 1];
        }
        dmvnorm_2d(x, y, mu, visigma, detsigma, n, dens);
        for (k = i; k <= N; k++)
        {
            sum[0] += dens[k - 1]*Xmat0[k - 1];
            sum[1] += dens[k - 1]*Xmat1[k - 1];
            sum[2] += dens[k - 1]*Xmat2[k - 1];
            sum[3] += dens[k - 1]*Xmat3[k - 1];
        }
    }

    free(x);
    free(y);
    free(dens);
    free(Xmat0);   
    free(Xmat1);
    free(Xmat2);    
    free(Xmat3);
}

/*************************************************************************
* Double sum used in density derivative estimation - 1st order
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* n - number of values
* sum - contains sum of density derivative values
*************************************************************************/

void dmvnormd1_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n,
    double *sum)
{
    int i, j, k, N;
    double sum1;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnormd1_2d(x, y, vsigma, r, n, derivt);

        for (k = i; k <= N; k++)
            sum1 = sum1 + derivt[k - 1];
    }

    sum[0] = sum1;
    free(x);
    free(y);
    free(derivt);
}


/*************************************************************************
* Double sum used in density derivative estimation - 2nd order
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* n - number of values
* sum - contains sum of density derivative values
*************************************************************************/

void dmvnormd2_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n,
    double *sum)
{
    int i, j, k, N;
    double sum1;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnormd2_2d(x, y, vsigma, r, n, derivt);

        for (k = i; k <= N; k++)
            sum1 = sum1 + derivt[k - 1];
    }

    sum[0] = sum1;
    free(x);
    free(y);
    free(derivt);
}



/*************************************************************************
* Double sum used in density derivative estimation - 3rd order
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* n - number of values
* sum - contains sum of density derivative values
*************************************************************************/

void dmvnormd3_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n,
    double *sum)
{
    int i, j, k, N;
    double sum1;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnormd3_2d(x, y, vsigma, r, n, derivt);  

        for (k = i; k <= N; k++)
            sum1 = sum1 + derivt[k - 1];
    }

    sum[0] = sum1;

    free(x);
    free(y);
    free(derivt);
}


/*************************************************************************
* Double sum used in density derivative estimation - 4th order
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* n - number of values
* sum - contains sum of density derivative values
*************************************************************************/

void dmvnormd4_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n,
    double *sum)
{
    int i, j, k, N;
    double sum1;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnormd4_2d(x, y, vsigma, r, n, derivt);

        for (k = i; k <= N; k++)
            sum1 = sum1 + derivt[k - 1];
    }

    sum[0] = sum1;
    free(x);
    free(y);
    free(derivt);
}

/*************************************************************************
* Double sum used in density derivative estimation - 5th order
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* n - number of values
* sum - contains sum of density derivative values
*************************************************************************/

void dmvnormd5_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n,
    double *sum)
{
    int i, j, k, N;
    double sum1;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnormd5_2d(x, y, vsigma, r, n, derivt);

        for (k = i; k <= N; k++)
            sum1 = sum1 + derivt[k - 1];
    }

    sum[0] = sum1;
    free(x);
    free(y);
    free(derivt);
}

/*************************************************************************
* Double sum used in density derivative estimation - 6th order
*
* Parameters
* x1 - x values
* x2 - y values
* v - vec Sigma
* r - (r1, r2) partial derivative
* n - number of values
* sum - contains sum of density derivative values
*************************************************************************/

void dmvnormd6_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n,
    double *sum)
{
    int i, j, k, N;
    double sum1;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    sum1 = 0.0;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
        }
        dmvnormd6_2d(x, y, vsigma, r, n, derivt);

        for (k = i; k <= N; k++)
            sum1 = sum1 + derivt[k - 1];
    }

    sum[0] = sum1;
    free(x);
    free(y);
    free(derivt);
}



/*************************************************************************
* Double sum of 4th order derivative normal density values * x x^T
*
* Parameters
* x1 - x values
* x2 - y values
* visigma - vec of inverse Sigma
* detsigma - det of Sigma
* n - number of values
* sum - contains sum of density values
*************************************************************************/

void dmvnormd4_2d_xxt_sum(double *x1, double *x2, double *visigma, 
			  int*r, int *n, double *sum)
{
    int i, j, k, N;
    double *Xmat0, *Xmat1, *Xmat2, *Xmat3;
    double *x, *y, *derivt;

    N = n[0];
    x = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    derivt = malloc(sizeof(double)*N);
    Xmat0 = malloc(sizeof(double)*N);
    Xmat1 = malloc(sizeof(double)*N);
    Xmat2 = malloc(sizeof(double)*N);
    Xmat3 = malloc(sizeof(double)*N);
     
    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            x[j - 1] = x1[i - 1] - x1[j - 1];
            y[j - 1] = x2[i - 1] - x2[j - 1];
            Xmat0[j - 1] = x[j - 1] * x[j - 1];
            Xmat2[j - 1] = Xmat1[j - 1] = x[j - 1] * y[j - 1];
            Xmat3[j - 1] = y[j - 1] * y[j - 1];
        }
        dmvnormd4_2d(x, y, visigma, r, n, derivt);
        for (k = i; k <= N; k++)
        {
            sum[0] += derivt[k - 1]*Xmat0[k - 1];
            sum[1] += derivt[k - 1]*Xmat1[k - 1];
            sum[2] += derivt[k - 1]*Xmat2[k - 1];
            sum[3] += derivt[k - 1]*Xmat3[k - 1];
        }
    }

    free(x);
    free(y);
    free(derivt);
    free(Xmat0);   
    free(Xmat1);
    free(Xmat2);    
    free(Xmat3);
}

/*************************************************************************
 * Bivariate t density
 *
 * Parameters
 * x - x values
 * y - y values
 * mu - mean
 * viSigma - vector form of inverse of dispersion matrix
 * 
 * n - number of values 
 * dens - contains density values
 *************************************************************************/
void dmvt_2d(double *x, double *y, double *mu, double *viSigma,  
	     double *df, int *n, double *dens)
{
  double xmu[2], d;
  int j;
  
  d = 2.0;
  for (j = 1; j <= n[0]; j++)
  {
    xmu[0] = x[j - 1] - mu[0];
    xmu[1] = y[j - 1] - mu[1];
    dens[j - 1] = pow( 1 + mult(xmu, viSigma, xmu, d)/df[0], -0.5*(-d +df[0]));
  } 
}

/*************************************************************************
* 3-dim normal density 
*************************************************************************/

void dmvnorm_3d(double *x1, double *x2, double *x3, double *mu, 
		double *visigma, double *detsigma, int *n, double *dens)
{
  double norm;
  double xmu[3];
  int j, d;
  
  d = 3;
  norm = 1/sqrt(pow(2*M_PI, d) * detsigma[0]);
  for (j = 1; j <= n[0]; j++)
  {
    xmu[0] = x1[j - 1] - mu[0];
    xmu[1] = x2[j - 1] - mu[1];
    xmu[2] = x3[j - 1] - mu[2]; 
    dens[j - 1] = norm * exp(-0.5 * mult(xmu, visigma, xmu, d));
  }
}

/*************************************************************************
* 3-dim normal density - double sum 
*************************************************************************/

void dmvnorm_3d_sum(double *x1, double *x2, double *x3, double *visigma, 
		    double *detsigma, int *n, double *sum)
{
  int i, j, k, N;
  double mu[] = {0.0, 0.0, 0.0};
  double sum1;
  double *y1, *y2, *y3, *dens;
  
  N = n[0];
  y1 = malloc(sizeof(double)*N);
  y2 = malloc(sizeof(double)*N);
  y3 = malloc(sizeof(double)*N);
  dens = malloc(sizeof(double)*N);
  sum1 = 0.0;
 
  for (i = 1; i <= N; i++)
  {
    for (j = i; j <= N; j++)
    {
      y1[j - 1] = x1[i - 1] - x1[j - 1];
      y2[j - 1] = x2[i - 1] - x2[j - 1];
      y3[j - 1] = x3[i - 1] - x3[j - 1];
    }
    dmvnorm_3d(y1, y2, y3, mu, visigma, detsigma, n, dens);
    for (k = i; k <= N; k++)
      sum1 = sum1 + dens[k - 1];
  }

  sum[0] = sum1;
  free(y1);
  free(y2);
  free(y3);
  free(dens);
}



/*************************************************************************
* 4-dim normal density 
*************************************************************************/

void dmvnorm_4d(double *x1, double *x2, double *x3, double *x4, double *mu, 
		double *visigma, double *detsigma, int *n, double *dens)
{
  double norm;
  double xmu[4];
  int j, d;
  
  d = 4;
  norm = 1/sqrt(pow(2*M_PI, d) * detsigma[0]);
  for (j = 1; j <= n[0]; j++)
  {
    xmu[0] = x1[j - 1] - mu[0];
    xmu[1] = x2[j - 1] - mu[1];
    xmu[2] = x3[j - 1] - mu[2]; 
    xmu[3] = x4[j - 1] - mu[3];
    dens[j - 1] = norm * exp(-0.5 * mult(xmu, visigma, xmu, d));
  }
}

/*************************************************************************
* 4-dim normal density - double sum 
*************************************************************************/

void dmvnorm_4d_sum(double *x1, double *x2, double *x3, double *x4, double *visigma, 
		    double *detsigma, int *n, double *sum)
{
  int i, j, k, N;
  double mu[] = {0.0, 0.0, 0.0, 0.0};
  double sum1;
  double *y1, *y2, *y3, *y4, *dens;
  
  N = n[0];
  y1 = malloc(sizeof(double)*N);
  y2 = malloc(sizeof(double)*N);
  y3 = malloc(sizeof(double)*N);
  y4 = malloc(sizeof(double)*N);
  dens = malloc(sizeof(double)*N);
  sum1 = 0.0;
 
  for (i = 1; i <= N; i++)
  {
    for (j = i; j <= N; j++)
    {
      y1[j - 1] = x1[i - 1] - x1[j - 1];
      y2[j - 1] = x2[i - 1] - x2[j - 1];
      y3[j - 1] = x3[i - 1] - x3[j - 1];
      y4[j - 1] = x4[i - 1] - x4[j - 1];
    }
    dmvnorm_4d(y1, y2, y3, y4, mu, visigma, detsigma, n, dens);
    for (k = i; k <= N; k++)
      sum1 = sum1 + dens[k - 1];
  }

  sum[0] = sum1;
  free(y1);
  free(y2);
  free(y3);
  free(y4);
  free(dens);
}


/*************************************************************************
* 2nd order partial derivative of 4-dim normal density 
*************************************************************************/

void dmvnormd2_4d(double *x1, double *x2, double *x3, double *x4, double *vsigma, 
		  int *r, int *n, double *derivt)
{
  double sigma11, sigma12, sigma13, sigma14, sigma22, sigma23, sigma24, 
    sigma33, sigma34, sigma44, y1, y2, y3, y4;
  int i;
    
  sigma11 = sqrt(vsigma[0]);
  sigma12 = vsigma[1];
  sigma13 = vsigma[2];
  sigma14 = vsigma[3];
  sigma22 = sqrt(vsigma[4]);
  sigma23 = vsigma[5];
  sigma24 = vsigma[6];
  sigma33 = sqrt(vsigma[7]);
  sigma34 = vsigma[8];
  sigma44 = sqrt(vsigma[9]);


  /* second order derivatives */
    
  for (i = 1; i <= n[0]; i++)
  {
    y1 = x1[i - 1];
    y2 = x2[i - 1];
    y3 = x3[i - 1];
    y4 = x4[i - 1];

  if ((r[0] == 2) && (r[1] == 0) && (r[2] == 0) && (r[3] == 0)) 	  
  {  
    derivt[i - 1] =
(-(1/(pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
              pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*
          pow(sigma11,2))) + pow(y1,2)/
      (pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
            pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*pow(sigma11,4)
       ))/(4.*pow(M_PI,2)*sigma11*sigma22*sigma33*sigma44);
  }
  else if ((r[0] == 0) && (r[1] == 2) && (r[2] == 0) && (r[3] == 0))
  {
    derivt[i - 1] =
      (-(1/(pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
              pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*
          pow(sigma22,2))) + pow(y2,2)/
      (pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
            pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*pow(sigma22,4)
       ))/(4.*pow(M_PI,2)*sigma11*sigma22*sigma33*sigma44);
  }
  else if ((r[0] == 0) && (r[1] == 0) && (r[2] == 2) && (r[3] == 0))
  {
    derivt[i - 1] = 
      (-(1/(pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
              pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*
          pow(sigma33,2))) + pow(y3,2)/
      (pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
            pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*pow(sigma33,4)
       ))/(4.*pow(M_PI,2)*sigma11*sigma22*sigma33*sigma44);
  }
  else if ((r[0] == 0) && (r[1] == 0) && (r[2] == 0) && (r[3] == 2))
  {    
    derivt[i - 1] =
	(-(1/(pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
              pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*
          pow(sigma44,2))) + pow(y4,2)/
      (pow(M_E,(pow(y1,2)/pow(sigma11,2) + pow(y2,2)/pow(sigma22,2) + 
            pow(y3,2)/pow(sigma33,2) + pow(y4,2)/pow(sigma44,2))/2.)*pow(sigma44,4)
       ))/(4.*pow(M_PI,2)*sigma11*sigma22*sigma33*sigma44);
  }
  else if ((r[0] == 1) && (r[1] == 1) && (r[2] == 0) && (r[3] == 0))
  {
    derivt[i - 1] = (y1*y2)/(4.*pow(M_E,(pow(y1,2)/pow(sigma11,2) + 
         pow(y2,2)/pow(sigma22,2) + pow(y3,2)/pow(sigma33,2) + 
         pow(y4,2)/pow(sigma44,2))/2.)*pow(M_PI,2)*pow(sigma11,3)*
			       pow(sigma22,3)*sigma33*sigma44);
  }
  else if ((r[0] == 1) && (r[1] == 0) && (r[2] == 1) && (r[3] == 0))
  {
    derivt[i - 1] = (y1*y3)/(4.*pow(M_E,(pow(y1,2)/pow(sigma11,2) + 
         pow(y2,2)/pow(sigma22,2) + pow(y3,2)/pow(sigma33,2) + 
         pow(y4,2)/pow(sigma44,2))/2.)*pow(M_PI,2)*pow(sigma11,3)*sigma22*
			     pow(sigma33,3)*sigma44);
  }
  else if ((r[0] == 1) && (r[1] == 0) && (r[2] == 0) && (r[3] == 1))
  {
    derivt[i - 1] = (y1*y4)/(4.*pow(M_E,(pow(y1,2)/pow(sigma11,2) + 
         pow(y2,2)/pow(sigma22,2) + pow(y3,2)/pow(sigma33,2) + 
         pow(y4,2)/pow(sigma44,2))/2.)*pow(M_PI,2)*pow(sigma11,3)*sigma22*
			     sigma33*pow(sigma44,3));
  }
  else if ((r[0] == 0) && (r[1] == 1) && (r[2] == 1) && (r[3] == 0))
  {    
    derivt[i - 1] = (y2*y3)/(4.*pow(M_E,(pow(y1,2)/pow(sigma11,2) + 
         pow(y2,2)/pow(sigma22,2) + pow(y3,2)/pow(sigma33,2) + 
         pow(y4,2)/pow(sigma44,2))/2.)*pow(M_PI,2)*sigma11*pow(sigma22,3)*
			     pow(sigma33,3)*sigma44);
  }
  else if ((r[0] == 0) && (r[1] == 1) && (r[2] == 0) && (r[3] == 1))		     
  {    
    derivt[i - 1] = (y2*y4)/(4.*pow(M_E,(pow(y1,2)/pow(sigma11,2) + 
         pow(y2,2)/pow(sigma22,2) + pow(y3,2)/pow(sigma33,2) + 
         pow(y4,2)/pow(sigma44,2))/2.)*pow(M_PI,2)*sigma11*pow(sigma22,3)*
			     sigma33*pow(sigma44,3));
  }
  else if ((r[0] == 0) && (r[1] == 0) && (r[2] == 1) && (r[3] == 1))
  {    
    derivt[i - 1] = (y3*y4)/(4.*pow(M_E,(pow(y1,2)/pow(sigma11,2) + 
         pow(y2,2)/pow(sigma22,2) + pow(y3,2)/pow(sigma33,2) + 
         pow(y4,2)/pow(sigma44,2))/2.)*pow(M_PI,2)*sigma11*sigma22*
			     pow(sigma33,3)*pow(sigma44,3));
  }
  }
}

/*************************************************************************
* 5-dim normal density 
*************************************************************************/

void dmvnorm_5d(double *x1, double *x2, double *x3, double *x4, double *x5, 
		double *mu, 
		double *visigma, double *detsigma, int *n, double *dens)
{
  double norm;
  double xmu[5];
  int j, d;
  
  d = 5;
  norm = 1/sqrt(pow(2*M_PI, d) * detsigma[0]);
  for (j = 1; j <= n[0]; j++)
  {
    xmu[0] = x1[j - 1] - mu[0];
    xmu[1] = x2[j - 1] - mu[1];
    xmu[2] = x3[j - 1] - mu[2]; 
    xmu[3] = x4[j - 1] - mu[3];
    xmu[4] = x5[j - 1] - mu[4]; 
    dens[j - 1] = norm * exp(-0.5 * mult(xmu, visigma, xmu, d));
  }
}

/*************************************************************************
* 5-dim normal density - double sum 
*************************************************************************/

void dmvnorm_5d_sum(double *x1, double *x2, double *x3, double *x4, double *x5, 
		    double *visigma, double *detsigma, int *n, double *sum)
{
  int i, j, k, N;
  double mu[] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double sum1;
  double *y1, *y2, *y3, *y4, *y5, *dens;
  
  N = n[0];
  y1 = malloc(sizeof(double)*N);
  y2 = malloc(sizeof(double)*N);
  y3 = malloc(sizeof(double)*N);
  y4 = malloc(sizeof(double)*N);
  y5 = malloc(sizeof(double)*N);
  dens = malloc(sizeof(double)*N);
  sum1 = 0.0;
 
  for (i = 1; i <= N; i++)
  {
    for (j = i; j <= N; j++)
    {
      y1[j - 1] = x1[i - 1] - x1[j - 1];
      y2[j - 1] = x2[i - 1] - x2[j - 1];
      y3[j - 1] = x3[i - 1] - x3[j - 1];
      y4[j - 1] = x4[i - 1] - x4[j - 1];
      y5[j - 1] = x5[i - 1] - x5[j - 1];
    }
    dmvnorm_5d(y1, y2, y3, y4, y5, mu, visigma, detsigma, n, dens);
    for (k = i; k <= N; k++)
      sum1 = sum1 + dens[k - 1];
  }

  sum[0] = sum1;
  free(y1);
  free(y2);
  free(y3);
  free(y4);
  free(y5);
  free(dens);
}





/*************************************************************************
* 6-dim normal density 
*************************************************************************/

void dmvnorm_6d(double *x1, double *x2, double *x3, double *x4, double *x5, 
		double *x6, double *mu, 
		double *visigma, double *detsigma, int *n, double *dens)
{
  double norm;
  double xmu[6];
  int j, d;
  
  d = 6;
  norm = 1/sqrt(pow(2*M_PI, d) * detsigma[0]);
  for (j = 1; j <= n[0]; j++)
  {
    xmu[0] = x1[j - 1] - mu[0];
    xmu[1] = x2[j - 1] - mu[1];
    xmu[2] = x3[j - 1] - mu[2]; 
    xmu[3] = x4[j - 1] - mu[3];
    xmu[4] = x5[j - 1] - mu[4]; 
    xmu[5] = x6[j - 1] - mu[5];
    dens[j - 1] = norm * exp(-0.5 * mult(xmu, visigma, xmu, d));
  }
}

/*************************************************************************
* 6-dim normal density - double sum 
*************************************************************************/

void dmvnorm_6d_sum(double *x1, double *x2, double *x3, double *x4, double *x5, 
		    double *x6, double *visigma, 
		    double *detsigma, int *n, double *sum)
{
  int i, j, k, N;
  double mu[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double sum1;
  double *y1, *y2, *y3, *y4, *y5, *y6, *dens;
  
  N = n[0];
  y1 = malloc(sizeof(double)*N);
  y2 = malloc(sizeof(double)*N);
  y3 = malloc(sizeof(double)*N);
  y4 = malloc(sizeof(double)*N);
  y5 = malloc(sizeof(double)*N);
  y6 = malloc(sizeof(double)*N);
  dens = malloc(sizeof(double)*N);
  sum1 = 0.0;
 
  for (i = 1; i <= N; i++)
  {
    for (j = i; j <= N; j++)
    {
      y1[j - 1] = x1[i - 1] - x1[j - 1];
      y2[j - 1] = x2[i - 1] - x2[j - 1];
      y3[j - 1] = x3[i - 1] - x3[j - 1];
      y4[j - 1] = x4[i - 1] - x4[j - 1];
      y5[j - 1] = x5[i - 1] - x5[j - 1];
      y6[j - 1] = x6[i - 1] - x6[j - 1];
    }
    dmvnorm_6d(y1, y2, y3, y4, y5, y6, mu, visigma, detsigma, n, dens);
    for (k = i; k <= N; k++)
      sum1 = sum1 + dens[k - 1];
  }

  sum[0] = sum1;
  free(y1);
  free(y2);
  free(y3);
  free(y4);
  free(y5);
  free(y6);
  free(dens);
}


