/* linbin3D.f -- translated by f2c (version 20050501).
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

/* ccccccccc FORTRAN subroutine linbin3D.f cccccccccc */
/* Obtains bin counts for trivariate data */
/* via the linear binning strategy. In this version */
/* observations outside the mesh are ignored. */
/* Last changed: 28 JUL 2005 */
/* Subroutine */ int lbthrd(double *x, int *n, double *a1, 
	double *a2, double *a3, double *b1, double *b2, 
	double *b3, int *m1, int *m2, int *m3, double *
	gcounts)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, li1, li2, li3, ind1, ind2, ind3, ind4, ind5, ind6, 
	    ind7, ind8;
    static double rem1, rem2, rem3, lxi1, lxi2, lxi3, delta1, delta2, 
	    delta3;

/*     Initialize grid counts to zero */
    /* Parameter adjustments */
    --gcounts;
    --x;

    /* Function Body */
    i__1 = *m1 * *m2 * *m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gcounts[i__] = 0.;
/* L10: */
    }
    delta1 = (*b1 - *a1) / (*m1 - 1);
    delta2 = (*b2 - *a2) / (*m2 - 1);
    delta3 = (*b3 - *a3) / (*m3 - 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lxi1 = (x[i__] - *a1) / delta1 + 1;
	lxi2 = (x[*n + i__] - *a2) / delta2 + 1;
	lxi3 = (x[(*n << 1) + i__] - *a3) / delta3 + 1;
/*        Find the integer part of "lxi1","lxi2" and "lxi3" */
	li1 = (int) lxi1;
	li2 = (int) lxi2;
	li3 = (int) lxi3;
	rem1 = lxi1 - li1;
	rem2 = lxi2 - li2;
	rem3 = lxi3 - li3;
	if (li1 >= 1) {
	    if (li1 < *m1) {
		if (li2 >= 1) {
		    if (li2 < *m2) {
			if (li3 >= 1) {
			    if (li3 < *m3) {
				ind1 = li1 + *m1 * (li2 - 1) + *m1 * *m2 * (
					li3 - 1);
				ind2 = li1 + 1 + *m1 * (li2 - 1) + *m1 * *m2 *
					 (li3 - 1);
				ind3 = li1 + *m1 * li2 + *m1 * *m2 * (li3 - 1)
					;
				ind4 = li1 + 1 + *m1 * li2 + *m1 * *m2 * (li3 
					- 1);
				ind5 = li1 + *m1 * (li2 - 1) + *m1 * *m2 * 
					li3;
				ind6 = li1 + 1 + *m1 * (li2 - 1) + *m1 * *m2 *
					 li3;
				ind7 = li1 + *m1 * li2 + *m1 * *m2 * li3;
				ind8 = li1 + 1 + *m1 * li2 + *m1 * *m2 * li3;
				gcounts[ind1] += (1 - rem1) * (1 - rem2) * (1 
					- rem3);
				gcounts[ind2] += rem1 * (1 - rem2) * (1 - 
					rem3);
				gcounts[ind3] += (1 - rem1) * rem2 * (1 - 
					rem3);
				gcounts[ind4] += rem1 * rem2 * (1 - rem3);
				gcounts[ind5] += (1 - rem1) * (1 - rem2) * 
					rem3;
				gcounts[ind6] += rem1 * (1 - rem2) * rem3;
				gcounts[ind7] += (1 - rem1) * rem2 * rem3;
				gcounts[ind8] += rem1 * rem2 * rem3;
			    }
			}
		    }
		}
	    }
	}
/* L20: */
    }
    return 0;
} /* lbthrd_ */

