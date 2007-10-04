/* linbin.f -- translated by f2c (version 20050501).
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

/* ccccccccc FORTRAN subroutine linbin.f cccccccccc */
/* Obtains bin counts for univariate data */
/* via the linear binning strategy. If "trun=0" then */
/* weight from end observations is given to corresponding */
/* end grid points. If "trun=1" then end observations */
/* are truncated. */
/* Last changed: 27/01/95 */
/* Subroutine */ int linbin(double *x, int *n, double *a, 
	double *b, int *m, int *trun, double *gcounts)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, li;
    static double rem, lxi, delta;

/*     Initialize grid counts to zero */
    /* Parameter adjustments */
    --gcounts;
    --x;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gcounts[i__] = 0.;
/* L10: */
    }
    delta = (*b - *a) / (*m - 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lxi = (x[i__] - *a) / delta + 1;
/*        Find integer part of "lxi" */
	li = (int) lxi;
	rem = lxi - li;
	if (li >= 1 && li < *m) {
	    gcounts[li] += 1 - rem;
	    gcounts[li + 1] += rem;
	} else if (li < 1 && *trun == 0) {
	    ++gcounts[1];
	} else if (li >= *m) {
	    if (li == *m || *trun == 0) {
		++gcounts[*m];
	    }
	}
/* L20: */
    }
    return 0;
} /* linbin_ */

