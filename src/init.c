#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void tseries_pp_sum (double* u, int* n, int* l, double* sum);
void tseries_quad_map (double *x, double *xi, double *a, int *n);
void tseries_arma (double *x, double *u, double *a, int *ar, int *ma,
		   int *arl, int *mal, int *max, int *n, int *intercept);
void tseries_bdstest_main (int *N, int *M, double *x, double *c,
			   double *cstan, double *EPS, int *TRACE);
void tseries_boot (double *x, double *xb, int *n, double *b, int *type);
void tseries_fit_garch (double *y, int *n, double *par, int *p, int *q,
			int *itmax, double *afctol, double *rfctol,
			double *xctol, double *xftol, double *fret,
			int *agrad, int *trace);
void tseries_ophess_garch (double *y, int *n, double *par, double *he,
			   int *p, int *q);
void tseries_pred_garch (double *y, double *h, int *n, double *par,
			 int *p, int *q, int *genuine);

static const R_CMethodDef CEntries[] = {
    {"tseries_pp_sum", (DL_FUNC) &tseries_pp_sum, 4},
    {"tseries_quad_map", (DL_FUNC) &tseries_quad_map, 4},
    {"tseries_arma", (DL_FUNC) &tseries_arma, 10},
    {"tseries_bdstest_main", (DL_FUNC) &tseries_bdstest_main, 7},
    {"tseries_boot", (DL_FUNC) &tseries_boot, 5},
    {"tseries_fit_garch", (DL_FUNC) &tseries_fit_garch, 13},
    {"tseries_ophess_garch", (DL_FUNC) &tseries_ophess_garch, 6},
    {"tseries_pred_garch", (DL_FUNC) &tseries_pred_garch, 7},
    {NULL, NULL, 0}
};

void R_init_tseries(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
