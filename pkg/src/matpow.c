/*
 *  Power of a matrix x^k := x x ... x, where x in an (n x n) matrix
 *  and k is an *integer*. Based on code originally written by Vincent
 *  Goulet for package actuar (and inspired from Octave in file
 *  .../src/xpow.cc) with slight shortcuts by Martin Maechler:
 */

#include "matpow.h"

/* .Call() this from R : */
SEXP R_matpow(SEXP x, SEXP k)
{
    if (!isMatrix(x))
	error(_("not a matrix"));

    SEXP dims = getAttrib(x, R_DimSymbol), z;
    int n = INTEGER(dims)[0],
	ktmp = INTEGER(k)[0]; /* need copy, as it is altered in matpow() */

/*     if (!isNumeric(x)) */
/* 	PROTECT(x = coerceVector(x, REALSXP)); /\* may give error... *\/ */
/*     else */
/* 	PROTECT(x = duplicate(x)); /\* since matpow() will alter it *\/ */

    PROTECT(x = coerceVector(x, REALSXP)); /* coercion to numeric */

    if (n != INTEGER(dims)[1]) {
	UNPROTECT(1);
	error(_("non-square matrix"));
    }
    if (n == 0)
	return(allocMatrix(REALSXP, 0, 0));

    PROTECT(z = allocMatrix(REALSXP, n, n));
    setAttrib(z, R_DimNamesSymbol,
	      getAttrib(x, R_DimNamesSymbol));

    matpow(REAL(x), n, ktmp, REAL(z));

    UNPROTECT(2);
    return z;
}

/* Compute z := x %^% k, x an (n x n) square "matrix" in column-order;
 * NB: x will be altered! The caller must make a copy if needed */
void matpow(double *x, int n, int k, double *z)
{
    if (k == 0) { /* return identity matrix */
	int i, j;
	for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
		z[i * n + j] = (i == j) ? 1.0 : 0.0;
	return;
    }
    else if (k < 0) {
	error(_("power must be a positive integer"));
    }
    else { /* k >= 1 */
	static const char *transa = "N";
	static const double one = 1.0, zero = 0.0;
	int nSqr = n * n;
	double /* temporary matrix */
	    *tmp  = (double *) R_alloc(nSqr, sizeof(double));

	/* Take powers in multiples of 2 until there is only one
	 * product left to make. That is, if k = 5, compute (x * x),
	 * then ((x * x) * (x * x)) and finally ((x * x) * (x * x)) *
	 * x. */
	Memcpy(z, x, (size_t) nSqr);

	k--;
	while (k > 0) {
	    if (k & 1) {	/* z := z * x */
		F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
				z, &n, x, &n, &zero, tmp, &n);
		Memcpy(z, tmp, (size_t) nSqr);
	    }
	    if(k == 1)
		break;
	    k >>= 1; /* efficient division by 2; now have k >= 1 */

	    /* x := x * x */
	    F77_CALL(dgemm)(transa, transa, &n, &n, &n, &one,
			    x, &n, x, &n, &zero, tmp, &n);
	    Memcpy(x, tmp, (size_t) nSqr);
	}
    }
}
