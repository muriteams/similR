// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// las
IntegerMatrix las(const ListOf< IntegerMatrix >& M, std::string rule, double threshold, bool self);
RcppExport SEXP _similR_las(SEXP MSEXP, SEXP ruleSEXP, SEXP thresholdSEXP, SEXP selfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ListOf< IntegerMatrix >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< std::string >::type rule(ruleSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type self(selfSEXP);
    rcpp_result_gen = Rcpp::wrap(las(M, rule, threshold, self));
    return rcpp_result_gen;
END_RCPP
}
// similarity
NumericMatrix similarity(const ListOf<IntegerMatrix>& M, const std::string& statistic, bool normalized);
RcppExport SEXP _similR_similarity(SEXP MSEXP, SEXP statisticSEXP, SEXP normalizedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ListOf<IntegerMatrix>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type statistic(statisticSEXP);
    Rcpp::traits::input_parameter< bool >::type normalized(normalizedSEXP);
    rcpp_result_gen = Rcpp::wrap(similarity(M, statistic, normalized));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_similR_las", (DL_FUNC) &_similR_las, 4},
    {"_similR_similarity", (DL_FUNC) &_similR_similarity, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_similR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
