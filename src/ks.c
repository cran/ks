#include <stdlib.h>
#include <math.h>
#define E  2.7182818284590452354

double mult(double *x, double *vecA, double *y, int d);
void mat_mult(double *vecA, double *vecB, int d, double *vecAB);

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
void dmvnormd4_2d_sum(double *x1, double *x2, double *vsigma, int *r, int *n, 
		      double *sum);
void dmvnormd4_2d_xxt_sum(double *x1, double *x2, double *visigma, 
			  int*r, int *n, double *sum);


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

void mat_mult(double *vecA, double *vecB, int d, double *vecAB)
{
  int i,j;
  int m = 0;

  for (j = 1; j <= d; j++)
    for (i = 1; i <= d; i++){
      m++;
      vecAB[m-1] = vecAB[m-1] + vecA[(i-1)*d + j]*vecB[(i-1)*d + j];
    }
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
	
      derivt[i - 1] = -(pow(M_E,(pow(sigma22,2)*pow(x,2) - 2*rho*sigma11*sigma22*x*y + 
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

