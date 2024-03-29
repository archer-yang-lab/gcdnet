/* GENERATED BY THE R FUNCTION CALL:
 * tools::package_native_routine_registration_skeleton("gcdnet", character_only = FALSE)
 */

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(erlassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hsvmlassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(loglassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lslassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sqsvmlassonet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"erlassonet",    (DL_FUNC) &F77_NAME(erlassonet),    26},
    {"hsvmlassonet",  (DL_FUNC) &F77_NAME(hsvmlassonet),  26},
    {"loglassonet",   (DL_FUNC) &F77_NAME(loglassonet),   25},
    {"lslassonet",    (DL_FUNC) &F77_NAME(lslassonet),    25},
    {"sqsvmlassonet", (DL_FUNC) &F77_NAME(sqsvmlassonet), 25},
    {NULL, NULL, 0}
};

void R_init_gcdnet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
