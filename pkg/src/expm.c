/*  ===== File part of R package expm =====
 *
 *  Function to compute the matrix exponential
 *
 *     exp(M) = sum(n = 0:Inf; M^n / n!),
 *
 *  where M is an (n x n) matrix.
 *
 *  The functions therein use LAPACK and BLAS routines. Nicely
 *  formatted man pages for these can be found at
 *
 *    <http://www.mathkeisan.com/UsersGuide/E/>
 *
 *  AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Christophe
 *  Dutang, based on code in package Matrix.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "locale.h"
#include "expm.h"


/* For matrix exponential calculations. Pade constants
 *
 *   n_{pqj} = [(p + q - j)! p!]/[(p + q)! j! (p - j)!]
 *
 * and
 *
 *   d_{pqj} = [(p + q - j)! q!]/[(p + q)! j! (q - j)!]
 *
 * for p = q = 8 and j = 1, ..., 8.
 */
const static double padec88 [] =
{
    5.0000000000000000e-1,
    1.1666666666666667e-1,
    1.6666666666666667e-2,
    1.6025641025641026e-3,
    1.0683760683760684e-4,
    4.8562548562548563e-6,
    1.3875013875013875e-7,
    1.9270852604185938e-9,
};


/* Matrix exponential exp(x), where x is an (n x n) matrix. Result z
 * is an (n x n) matrix. Mostly lifted from the core of fonction
 * expm() of package Matrix, which is itself based on the function of
 * the same name in Octave.
 */
void expm(double *x, int n, double *z)
{
    if (n == 1)
	z[0] = exp(x[0]);		/* scalar exponential */
    else
    {
	/* Constants */
	int i, j;
	int nsqr = n * n, np1 = n + 1, is_uppertri = TRUE;
	int iloperm, ihiperm, iloscal, ihiscal, info, sqrpowscal;
	int lwork, is_conjug, is_diag;
	double infnorm, trshift, one = 1.0, zero = 0.0, m1pj = -1;
	double tmp;
	char jobVL[1], jobVR[1];
	char *transa = "N";
	Rcomplex cone, czero;

	/* Arrays */
	int *pivot    = (int *) R_alloc(n, sizeof(int)); /* pivot vector */
	int *invperm  = (int *) R_alloc(n, sizeof(int)); /* inverse permutation vector */
	double *perm  = (double *) R_alloc(n, sizeof(double)); /* permutation array */
	double *scale = (double *) R_alloc(n, sizeof(double)); /* scale array */
	double *work  = (double *) R_alloc(nsqr, sizeof(double)); /* workspace array */
	double *npp   = (double *) R_alloc(nsqr, sizeof(double)); /* num. power Pade */
	double *dpp   = (double *) R_alloc(nsqr, sizeof(double)); /* denom. power Pade */
	double *left, *right, *workdiag; /* left and right eigenvectors and workspace for diagonalisation */
	double *wR = (double *) R_alloc(n, sizeof(double)); /* real part of eigenvalues */
    double *wI = (double *) R_alloc(n, sizeof(double)); /* imaginary part of eigenvalues */
	Rcomplex *eigvect = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* (right) eigenvectors matrix */
	Rcomplex *eigvectinv = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* its inverse */
	Rcomplex *expeigval; /* complex matrix diag(exp(eigenvalues)) */
	Rcomplex *ctmp;
	int *ipiv = (int *) R_alloc(n, sizeof(int)); /* permutation vector */
	R_CheckStack();

	Memcpy(z, x, nsqr);

	/* Check if matrix x is upper triangular; stop checking as
	 * soon as a non-zero value is found below the diagonal. */
	for (i = 0; i < n - 1 && is_uppertri; i++)
	    for (j = i + 1; j < n; j++)
		if (!(is_uppertri = x[i * n + j] == 0.0))
		    break;
	
	/* Test if x is diagonalisable by computing its eigenvalues and (right) eigenvectors */
	/* base on modla_rg in lapack.c, used in eigen.R */
	jobVL[0] = 'N';
    left = (double *) 0;
	jobVR[0] = 'V';
	right = (double *) R_alloc(nsqr, sizeof(double));

	/* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgeev)(jobVL, jobVR, &n, z, &n, wR, wI, left, &n, right, &n, &tmp, &lwork, &info);
	if (info != 0)
		error(_("error code %d from Lapack routine dgeev"), info); 
    lwork = (int) tmp;
    workdiag = (double *) R_alloc(lwork, sizeof(double)); 

	/* compute eigenvalues and (right) eigenvectors */
	F77_CALL(dgeev)(jobVL, jobVR, &n, z, &n, wR, wI, left, &n, right, &n, workdiag, &lwork, &info); 
    if (info != 0)
		error(_("error code %d from Lapack routine dgeev"), info);
			
	/* try to invert the eigenvectors matrix */	
	/* build the Rcomplex matrix with eigenvectors */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{					
			is_conjug = 0;
			if(i < n-1)
			{	/* conjugate eigenvalues */
				if(wR[i] == wR[i+1] && wI[i] == -wI[i+1])
				{
					is_conjug = 1;
					eigvect[i * n + j].r = right[i * n + j];
					eigvect[i * n + j].i = right[(i+1) * n + j];
				}	
			}
			if(i > 0)
			{	/* conjugate eigenvalues */
				if(wR[i] == wR[i-1] && wI[i] == -wI[i-1])
				{
					is_conjug = 1;
					eigvect[i * n + j].r = right[(i-1) * n + j];
					eigvect[i * n + j].i = -right[i * n + j];					
				}	
			}
			/* real eigenvalues */
			if(!is_conjug)
			{
				eigvect[i * n + j].r = right[i * n + j];
				eigvect[i * n + j].i = 0.0;
			}
			/* eigvectinv initialise with the identity matrix */
			eigvectinv[i * n +j].r = (i == j) ? 1.0 : 0.0;
			eigvectinv[i * n +j].i = 0.0;
		}
	}
		
	/* solve a linear complex equation system with eigvectinv equals 
	 * to matrix identity. hence, on exit eigvectinv contains the 
	 * inverse of complex matrix eigvect */
	/* base on solve.R */ 
	F77_CALL(zgesv)(&n, &n, eigvect, &n, ipiv, eigvectinv, &n, &info); 
	
	
    if (info < 0)
	error(_("argument %d of Lapack routine dgesv had invalid value"), -info);
	if (info == 0) /* x is diagonalisable */
	{
		is_diag = 1;
		/* expeigval is the complex matrix diag(exp(eigenvalues)) */
		/* base on z_exp in complex.c */
		expeigval = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex));
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if(i == j)
				{
					expeigval[i * n +j].r = exp(wR[i]) * cos(wI[i]);
					expeigval[i * n +j].i = exp(wR[i]) * sin(wI[i]);
				}
				else
				{
					expeigval[i * n +j].r = 0.0;
					expeigval[i * n +j].i = 0.0;
				}
			}	
		}
		
	
		cone.r = 1.0; cone.i = czero.r = czero.i = 0.0;
		ctmp = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex));
		/* compute (complex) matrix product: tmp <- eigvect * expeigval */
		F77_CALL(zgemm)(transa, transa, &n, &n, &n, &cone, eigvect, &n, expeigval, &n, &czero, ctmp, &n);
		
		/* compute (complex) matrix product: expeigval <- tmp * eigvectinv */
		F77_CALL(zgemm)(transa, transa, &n, &n, &n, &cone, ctmp, &n, eigvectinv, &n, &czero, expeigval, &n);
		
		/* store the real part in z (NB the matrix exponential is real) */
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) 
				z[i * n + j] = expeigval[i * n + j].r;
		
	}
	else
		is_diag = 0;				
			
	if(!is_diag)
	{
	
		/* Step 1 of preconditioning: shift diagonal by average
		 * diagonal if positive. */
		trshift = 0.0;
		for (i = 0; i < n; i++)
			trshift += x[i * np1];
		trshift /= n;		/* average diagonal element */
		if (trshift > 0.0)
			for (i = 0; i < n; i++)
			z[i * np1] -= trshift;

		/* Step 2 of preconditioning: balancing with dgebal. */
		if (is_uppertri)
		{
			/* no need to permute if x is upper triangular */
			iloperm = 1;
			ihiperm = n;
		}
		else
		{
			F77_CALL(dgebal)("P", &n, z, &n, &iloperm, &ihiperm, perm, &info);
			if (info)
			error(_("LAPACK routine dgebal returned info code %d when permuting"), info);
		}
		F77_CALL(dgebal)("S", &n, z, &n, &iloscal, &ihiscal, scale, &info);
		if (info)
			error(_("LAPACK routine dgebal returned info code %d when scaling"), info);

		/* Step 3 of preconditioning: Scaling according to infinity
		 * norm (a priori always needed). */
		infnorm = F77_CALL(dlange)("I", &n, &n, z, &n, work);
		sqrpowscal = (infnorm > 0) ? imax2((int) 1 + log(infnorm)/log(2.0), 0) : 0;
		if (sqrpowscal > 0)
		{
			double scalefactor = R_pow_di(2, sqrpowscal);
			for (i = 0; i < nsqr; i++)
			z[i] /= scalefactor;
		}

		/* Pade approximation (p = q = 8): compute x^8, x^7, x^6,
		 * ..., x^1 */
		for (i = 0; i < nsqr; i++)
		{
			npp[i] = 0.0;
			dpp[i] = 0.0;
		}
		for (j = 7; j >= 0; j--)
		{
			/* npp = z * npp + padec88[j] * z */
			F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, npp,
					 &n, &zero, work, &n);
			/* npp <- work + padec88[j] * z */
			for (i = 0; i < nsqr; i++)
			npp[i] = work[i] + padec88[j] * z[i];
			/* dpp = z * dpp + (-1)^j * padec88[j] * z */
			F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, dpp,
					 &n, &zero, work, &n);
			for (i = 0; i < nsqr; i++)
			dpp[i] = work[i] + m1pj * padec88[j] * z[i];
			m1pj *= -1;		/* (-1)^j */
		}
		/* power 0 */
		for (i = 0; i < nsqr; i++)
			dpp[i] *= -1.0;
		for (j = 0; j < n; j++)
		{
			npp[j * np1] += 1.0;
			dpp[j * np1] += 1.0;
		}

		/* Pade approximation is (dpp)^-1 * npp. */
		F77_CALL(dgetrf) (&n, &n, dpp, &n, pivot, &info);
		if (info)
			error(_("LAPACK routine dgetrf returned info code %d"), info);
		F77_CALL(dgetrs) ("N", &n, &n, dpp, &n, pivot, npp, &n, &info);
		if (info)
			error(_("LAPACK routine dgetrs returned info code %d"), info);

		Memcpy(z, npp, nsqr);

		/* Now undo all of the preconditioning */
		/* Preconditioning 3: square the result for every power of 2 */
		while (sqrpowscal--)
		{
			F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, z, &n,
					z, &n, &zero, work, &n);
			Memcpy(z, work, nsqr);
		}
		/* Preconditioning 2: apply inverse scaling */
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			z[i + j * n] *= scale[i]/scale[j];

		/* Inverse permuation if x is not upper triangular and 'perm'
		 * is not the identity permutation */
		if ((iloperm != 1 || ihiperm != n) && !is_uppertri)
		{
			/* balancing permutation vector */
			for (i = 0; i < n; i++)
			invperm[i] = i;	/* identity permutation */

			/* leading permutations applied in forward order */
			for (i = 0; i < (iloperm - 1); i++)
			{
			int permutedindex = (int) (perm[i]) - 1;
			int tmp = invperm[i];
			invperm[i] = invperm[permutedindex];
			invperm[permutedindex] = tmp;
			}

			/* trailing permutations applied in reverse order */
			for (i = n - 1; i >= ihiperm; i--)
			{
			int permutedindex = (int) (perm[i]) - 1;
			int tmp = invperm[i];
			invperm[i] = invperm[permutedindex];
			invperm[permutedindex] = tmp;
			}

			/* construct inverse balancing permutation vector */
			Memcpy(pivot, invperm, n);
			for (i = 0; i < n; i++)
			invperm[pivot[i]] = i;

			/* apply inverse permutation */
			Memcpy(work, z, nsqr);
			for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				z[i + j * n] = work[invperm[i] + invperm[j] * n];
		}
		/* Preconditioning 1: Trace normalization */
		if (trshift > 0)
		{
			double mult = exp(trshift);
			for (i = 0; i < nsqr; i++)
			z[i] *= mult;
		}
	}
    }
}

/* Main function, the only one used by .External(). */
SEXP do_expm(SEXP x)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;

    if (!isNumeric(x) || !isMatrix(x))
	error(_("invalid argument"));

    dims = getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m)
	error(_("non-square matrix"));
    if (n == 0)
	return(allocVector(REALSXP, 0));

    PROTECT(z = allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    expm(rx, n, rz);

    setAttrib(z, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(1);

    return z;
}
