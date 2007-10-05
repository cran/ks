/* linbin2D.f -- translated by f2c (version 20050501).
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

/* ccccccccc FORTRAN subroutine linbin2D.f cccccccccc */
/* Obtains bin counts for bivariate data */
/* via the linear binning strategy. In this version */
/* observations outside the mesh are ignored. */
/* Last changed: 25 AUG 1995 */
/* Subroutine */ int lbtwod(double *x, int  *n, double *a1, 
	double *a2, double *b1, double *b2, int *m1, int *
	m2, double *gcounts)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, li1, li2, ind1, ind2, ind3, ind4;
    static double rem1, rem2, lxi1, lxi2, delta1, delta2;

/*     Initialize grid counts to zero */
    /* Parameter adjustments */
    --gcounts;
    --x;

    /* Function Body */
    i__1 = *m1 * *m2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gcounts[i__] = 0.;
/* L10: */
    }
    delta1 = (*b1 - *a1) / (*m1 - 1);
    delta2 = (*b2 - *a2) / (*m2 - 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lxi1 = (x[i__] - *a1) / delta1 + 1;
	lxi2 = (x[*n + i__] - *a2) / delta2 + 1;
/*        Find the integer part of "lxi1" and "lxi2" */
	li1 = (int) lxi1;
	li2 = (int) lxi2;
	rem1 = lxi1 - li1;
	rem2 = lxi2 - li2;
	if (li1 >= 1) {
	    if (li2 >= 1) {
		if (li1 < *m1) {
		    if (li2 < *m2) {
			ind1 = *m1 * (li2 - 1) + li1;
			ind2 = *m1 * li2 + li1;
			ind3 = *m1 * (li2 - 1) + li1 + 1;
			ind4 = *m1 * li2 + li1 + 1;
			gcounts[ind1] += (1 - rem1) * (1 - rem2);
			gcounts[ind2] += rem1 * (1 - rem2);
			gcounts[ind3] += (1 - rem1) * rem2;
			gcounts[ind4] += rem1 * rem2;
		    }
		}
	    }
	}
/* L20: */
    }
    return 0;
} /* lbtwod_ */

