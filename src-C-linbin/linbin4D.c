/* linbin4D.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* ccccccccc FORTRAN subroutine linbin4D.f cccccccccc */
/* Obtains bin counts for quadrivariate data */
/* via the linear binning strategy. In this version */
/* observations outside the mesh are ignored. */
/* Last changed: 31 AUG 2005 */
/* Subroutine */ int lbfoud(double *x, int *n, double *a1, 
	double *a2, double *a3, double *a4, double *b1, 
	double *b2, double *b3, double *b4, int *m1, int *
	m2, int *m3, int *m4, double *gcounts)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, li1, li2, li3, li4, ind1, ind2, ind3, ind4, ind5, 
	    ind6, ind7, ind8, ind9;
    static double rem1, rem2, rem3, rem4, lxi1, lxi2, lxi3, lxi4;
    static int ind10, ind11, ind12, ind13, ind14, ind15, ind16;
    static double delta1, delta2, delta3, delta4;

/*     Initialize grid counts to zero */
    /* Parameter adjustments */
    --gcounts;
    --x;

    /* Function Body */
    i__1 = *m1 * *m2 * *m3 * *m4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gcounts[i__] = 0.;
/* L10: */
    }
    delta1 = (*b1 - *a1) / (*m1 - 1);
    delta2 = (*b2 - *a2) / (*m2 - 1);
    delta3 = (*b3 - *a3) / (*m3 - 1);
    delta4 = (*b4 - *a4) / (*m4 - 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lxi1 = (x[i__] - *a1) / delta1 + 1;
	lxi2 = (x[*n + i__] - *a2) / delta2 + 1;
	lxi3 = (x[(*n << 1) + i__] - *a3) / delta3 + 1;
	lxi4 = (x[*n * 3 + i__] - *a4) / delta4 + 1;
/*        Find the integer part of "lxi1","lxi2", "lxi3" and "lxi4" */
	li1 = (int) lxi1;
	li2 = (int) lxi2;
	li3 = (int) lxi3;
	li4 = (int) lxi4;
	rem1 = lxi1 - li1;
	rem2 = lxi2 - li2;
	rem3 = lxi3 - li3;
	rem4 = lxi4 - li4;
	if (li1 >= 1) {
	    if (li1 < *m1) {
		if (li2 >= 1) {
		    if (li2 < *m2) {
			if (li3 >= 1) {
			    if (li3 < *m3) {
				if (li4 > 1) {
				    if (li4 < *m4) {
					ind1 = li1 + *m1 * (li2 - 1) + *m1 * *
						m2 * (li3 - 1) + *m1 * *m2 * *
						m3 * (li4 - 1);
					ind2 = li1 + 1 + *m1 * (li2 - 1) + *
						m1 * *m2 * (li3 - 1) + *m1 * *
						m2 * *m3 * (li4 - 1);
					ind3 = li1 + *m1 * li2 + *m1 * *m2 * (
						li3 - 1) + *m1 * *m2 * *m3 * (
						li4 - 1);
					ind4 = li1 + 1 + *m1 * li2 + *m1 * *
						m2 * (li3 - 1) + *m1 * *m2 * *
						m3 * (li4 - 1);
					ind5 = li1 + *m1 * (li2 - 1) + *m1 * *
						m2 * li3 + *m1 * *m2 * *m3 * (
						li4 - 1);
					ind6 = li1 + 1 + *m1 * (li2 - 1) + *
						m1 * *m2 * li3 + *m1 * *m2 * *
						m3 * (li4 - 1);
					ind7 = li1 + *m1 * li2 + *m1 * *m2 * 
						li3 + *m1 * *m2 * *m3 * (li4 
						- 1);
					ind8 = li1 + 1 + *m1 * li2 + *m1 * *
						m2 * li3 + *m1 * *m2 * *m3 * (
						li4 - 1);
					ind9 = li1 + *m1 * (li2 - 1) + *m1 * *
						m2 * (li3 - 1) + *m1 * *m2 * *
						m3 * li4;
					ind10 = li1 + 1 + *m1 * (li2 - 1) + *
						m1 * *m2 * (li3 - 1) + *m1 * *
						m2 * *m3 * li4;
					ind11 = li1 + *m1 * li2 + *m1 * *m2 * 
						(li3 - 1) + *m1 * *m2 * *m3 * 
						li4;
					ind12 = li1 + 1 + *m1 * li2 + *m1 * *
						m2 * (li3 - 1) + *m1 * *m2 * *
						m3 * li4;
					ind13 = li1 + *m1 * (li2 - 1) + *m1 * 
						*m2 * li3 + *m1 * *m2 * *m3 * 
						li4;
					ind14 = li1 + 1 + *m1 * (li2 - 1) + *
						m1 * *m2 * li3 + *m1 * *m2 * *
						m3 * li4;
					ind15 = li1 + *m1 * li2 + *m1 * *m2 * 
						li3 + *m1 * *m2 * *m3 * li4;
					ind16 = li1 + 1 + *m1 * li2 + *m1 * *
						m2 * li3 + *m1 * *m2 * *m3 * 
						li4;
					gcounts[ind1] += (1 - rem1) * (1 - 
						rem2) * (1 - rem3) * (1 - 
						rem4);
					gcounts[ind2] += rem1 * (1 - rem2) * (
						1 - rem3) * (1 - rem4);
					gcounts[ind3] += (1 - rem1) * rem2 * (
						1 - rem3) * (1 - rem4);
					gcounts[ind4] += rem1 * rem2 * (1 - 
						rem3) * (1 - rem4);
					gcounts[ind5] += (1 - rem1) * (1 - 
						rem2) * rem3 * (1 - rem4);
					gcounts[ind6] += rem1 * (1 - rem2) * 
						rem3 * (1 - rem4);
					gcounts[ind7] += (1 - rem1) * rem2 * 
						rem3 * (1 - rem4);
					gcounts[ind8] += rem1 * rem2 * rem3 * 
						(1 - rem4);
					gcounts[ind9] += (1 - rem1) * (1 - 
						rem2) * (1 - rem3) * rem4;
					gcounts[ind10] += rem1 * (1 - rem2) * 
						(1 - rem3) * rem4;
					gcounts[ind11] += (1 - rem1) * rem2 * 
						(1 - rem3) * rem4;
					gcounts[ind12] += rem1 * rem2 * (1 - 
						rem3) * rem4;
					gcounts[ind13] += (1 - rem1) * (1 - 
						rem2) * rem3 * rem4;
					gcounts[ind14] += rem1 * (1 - rem2) * 
						rem3 * rem4;
					gcounts[ind15] += (1 - rem1) * rem2 * 
						rem3 * rem4;
					gcounts[ind16] += rem1 * rem2 * rem3 *
						 rem4;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
/* L20: */
    }
    return 0;
} /* lbfoud_ */

