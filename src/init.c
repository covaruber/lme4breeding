#include "pedigree.h"
#include <R_ext/Rdynload.h>
static R_CallMethodDef CallEntries[] = {
    {"pedigree_chol", (DL_FUNC) &pedigree_chol, 2},
    {"pedigree_inbreeding", (DL_FUNC) &pedigree_inbreeding, 1},
    {NULL, NULL, 0}
};

/** Initializer for the pedigreemm package.
 *
 *  Register routines that can be called directly from R.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_pedigreemm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

