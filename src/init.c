#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP scoringRules_auxcrpsC(SEXP, SEXP);
extern SEXP scoringRules_bvarFcstC(SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_carterkohn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_crpsmixnC(SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_dmixnC(SEXP, SEXP, SEXP);
extern SEXP scoringRules_dnormC(SEXP);
extern SEXP scoringRules_drawbetaC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_drawMultinomC(SEXP);
extern SEXP scoringRules_drawsigmaC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_energyscoreC(SEXP, SEXP);
extern SEXP scoringRules_euclnormC(SEXP);
extern SEXP scoringRules_filterMarkovMixtureC(SEXP, SEXP, SEXP);
extern SEXP scoringRules_getfcsts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_lsmixnC(SEXP, SEXP, SEXP, SEXP);
extern SEXP scoringRules_makeregs_fcC(SEXP, SEXP);
extern SEXP scoringRules_matmult(SEXP, SEXP);
extern SEXP scoringRules_meye(SEXP);
extern SEXP scoringRules_mvndrawC(SEXP, SEXP);
extern SEXP scoringRules_pmixnC(SEXP, SEXP, SEXP);
extern SEXP scoringRules_pnormC(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"scoringRules_auxcrpsC",             (DL_FUNC) &scoringRules_auxcrpsC,             2},
    {"scoringRules_bvarFcstC",            (DL_FUNC) &scoringRules_bvarFcstC,            4},
    {"scoringRules_carterkohn",           (DL_FUNC) &scoringRules_carterkohn,           9},
    {"scoringRules_crpsmixnC",            (DL_FUNC) &scoringRules_crpsmixnC,            4},
    {"scoringRules_dmixnC",               (DL_FUNC) &scoringRules_dmixnC,               3},
    {"scoringRules_dnormC",               (DL_FUNC) &scoringRules_dnormC,               1},
    {"scoringRules_drawbetaC",            (DL_FUNC) &scoringRules_drawbetaC,            5},
    {"scoringRules_drawMultinomC",        (DL_FUNC) &scoringRules_drawMultinomC,        1},
    {"scoringRules_drawsigmaC",           (DL_FUNC) &scoringRules_drawsigmaC,           9},
    {"scoringRules_energyscoreC",         (DL_FUNC) &scoringRules_energyscoreC,         2},
    {"scoringRules_euclnormC",            (DL_FUNC) &scoringRules_euclnormC,            1},
    {"scoringRules_filterMarkovMixtureC", (DL_FUNC) &scoringRules_filterMarkovMixtureC, 3},
    {"scoringRules_getfcsts",             (DL_FUNC) &scoringRules_getfcsts,             6},
    {"scoringRules_lsmixnC",              (DL_FUNC) &scoringRules_lsmixnC,              4},
    {"scoringRules_makeregs_fcC",         (DL_FUNC) &scoringRules_makeregs_fcC,         2},
    {"scoringRules_matmult",              (DL_FUNC) &scoringRules_matmult,              2},
    {"scoringRules_meye",                 (DL_FUNC) &scoringRules_meye,                 1},
    {"scoringRules_mvndrawC",             (DL_FUNC) &scoringRules_mvndrawC,             2},
    {"scoringRules_pmixnC",               (DL_FUNC) &scoringRules_pmixnC,               3},
    {"scoringRules_pnormC",               (DL_FUNC) &scoringRules_pnormC,               1},
    {NULL, NULL, 0}
};

void R_init_scoringRules(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
