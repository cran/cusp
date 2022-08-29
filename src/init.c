#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .External calls */
extern SEXP cuspnc(SEXP);

static const R_ExternalMethodDef ExternalEntries[] = {
  {"cuspnc", (DL_FUNC) &cuspnc, 7},
  {NULL, NULL, 0}
};

void R_init_cusp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
  R_useDynamicSymbols(dll, FALSE);
}

