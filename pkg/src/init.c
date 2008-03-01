/*
 *  Native routines registration
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "expm-eigen.h"
#include "expm.h"
#include "matpow.h"

static const R_CallMethodDef CallEntries[] = {
    {"do_expm", (DL_FUNC) &do_expm, 2},
    {"R_matpow", (DL_FUNC) &R_matpow, 2},
    {"do_expm_eigen", (DL_FUNC) &do_expm_eigen, 2},
    {NULL, NULL, 0}
};

void R_init_expm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    /* callable C code from other packages C code :*/
    R_RegisterCCallable("expm", "expm", (DL_FUNC) expm);
    R_RegisterCCallable("matpow", "matpow", (DL_FUNC) matpow);
    R_RegisterCCallable("expm_eigen", "expm_eigen", (DL_FUNC) expm_eigen);
}
