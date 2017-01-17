#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_pp_sum (double* u, int* n, int* l, double* sum);
void R_quad_map (double *x, double *xi, double *a, int *n);
void R_arma (double *x, double *u, double *a, int *ar, int *ma, int *arl, 
	     int *mal, int *max, int *n, int *intercept);
void R_bdstest_main (int *N, int *M, double *x, double *c, double *cstan,
		     double *EPS, int *TRACE);
void R_boot (double *x, double *xb, int *n, double *b, int *type);
void R_fit_garch (double *y, int *n, double *par, int *p, int *q, int *itmax, 
		  double *afctol, double *rfctol, double *xctol, double *xftol,
		  double *fret, int *agrad, int *trace);
void R_ophess_garch (double *y, int *n, double *par, double *he, int *p, int *q);
void R_pred_garch (double *y, double *h, int *n, double *par, int *p, int *q, int *genuine);

static const R_CMethodDef CEntries[] = {
    {"R_pp_sum", (DL_FUNC) &R_pp_sum, 4},
    {"R_quad_map", (DL_FUNC) &R_quad_map, 4},
    {"R_arma", (DL_FUNC) &R_arma, 10},
    {"R_bdstest_main", (DL_FUNC) &R_bdstest_main, 7},
    {"R_boot", (DL_FUNC) &R_boot, 5},
    {"R_fit_garch", (DL_FUNC) &R_fit_garch, 13},
    {"R_ophess_garch", (DL_FUNC) &R_ophess_garch, 6},
    {"R_pred_garch", (DL_FUNC) &R_pred_garch, 7},
    {NULL, NULL, 0}
};

void R_init_tseries(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
