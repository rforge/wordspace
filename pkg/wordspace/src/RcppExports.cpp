// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CPP_col_dist_dense
NumericVector CPP_col_dist_dense(NumericMatrix x, NumericMatrix y, int metric_code, double param1, bool symmetric);
RcppExport SEXP wordspace_CPP_col_dist_dense(SEXP xSEXP, SEXP ySEXP, SEXP metric_codeSEXP, SEXP param1SEXP, SEXP symmetricSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP );
        Rcpp::traits::input_parameter< int >::type metric_code(metric_codeSEXP );
        Rcpp::traits::input_parameter< double >::type param1(param1SEXP );
        Rcpp::traits::input_parameter< bool >::type symmetric(symmetricSEXP );
        NumericVector __result = CPP_col_dist_dense(x, y, metric_code, param1, symmetric);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CPP_row_norms_dense
NumericVector CPP_row_norms_dense(NumericMatrix x, int norm_code, double p_norm = 2.0);
RcppExport SEXP wordspace_CPP_row_norms_dense(SEXP xSEXP, SEXP norm_codeSEXP, SEXP p_normSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< int >::type norm_code(norm_codeSEXP );
        Rcpp::traits::input_parameter< double >::type p_norm(p_normSEXP );
        NumericVector __result = CPP_row_norms_dense(x, norm_code, p_norm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CPP_row_norms_sparse
NumericVector CPP_row_norms_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector x, int norm_code, double p_norm = 2.0);
RcppExport SEXP wordspace_CPP_row_norms_sparse(SEXP nrSEXP, SEXP ncSEXP, SEXP pSEXP, SEXP row_ofSEXP, SEXP xSEXP, SEXP norm_codeSEXP, SEXP p_normSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nr(nrSEXP );
        Rcpp::traits::input_parameter< int >::type nc(ncSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type row_of(row_ofSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< int >::type norm_code(norm_codeSEXP );
        Rcpp::traits::input_parameter< double >::type p_norm(p_normSEXP );
        NumericVector __result = CPP_row_norms_sparse(nr, nc, p, row_of, x, norm_code, p_norm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CPP_col_norms_dense
NumericVector CPP_col_norms_dense(NumericMatrix x, int norm_code, double p_norm = 2.0);
RcppExport SEXP wordspace_CPP_col_norms_dense(SEXP xSEXP, SEXP norm_codeSEXP, SEXP p_normSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< int >::type norm_code(norm_codeSEXP );
        Rcpp::traits::input_parameter< double >::type p_norm(p_normSEXP );
        NumericVector __result = CPP_col_norms_dense(x, norm_code, p_norm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CPP_col_norms_sparse
NumericVector CPP_col_norms_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector x, int norm_code, double p_norm = 2.0);
RcppExport SEXP wordspace_CPP_col_norms_sparse(SEXP nrSEXP, SEXP ncSEXP, SEXP pSEXP, SEXP row_ofSEXP, SEXP xSEXP, SEXP norm_codeSEXP, SEXP p_normSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nr(nrSEXP );
        Rcpp::traits::input_parameter< int >::type nc(ncSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type row_of(row_ofSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< int >::type norm_code(norm_codeSEXP );
        Rcpp::traits::input_parameter< double >::type p_norm(p_normSEXP );
        NumericVector __result = CPP_col_norms_sparse(nr, nc, p, row_of, x, norm_code, p_norm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CPP_dsm_score_dense
NumericMatrix CPP_dsm_score_dense(NumericMatrix f, NumericVector f1, NumericVector f2, double N, int am_code, int sparse, int transform_code);
RcppExport SEXP wordspace_CPP_dsm_score_dense(SEXP fSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP NSEXP, SEXP am_codeSEXP, SEXP sparseSEXP, SEXP transform_codeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type f(fSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type f1(f1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type f2(f2SEXP );
        Rcpp::traits::input_parameter< double >::type N(NSEXP );
        Rcpp::traits::input_parameter< int >::type am_code(am_codeSEXP );
        Rcpp::traits::input_parameter< int >::type sparse(sparseSEXP );
        Rcpp::traits::input_parameter< int >::type transform_code(transform_codeSEXP );
        NumericMatrix __result = CPP_dsm_score_dense(f, f1, f2, N, am_code, sparse, transform_code);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CPP_dsm_score_sparse
NumericVector CPP_dsm_score_sparse(int nr, int nc, IntegerVector p, IntegerVector row_of, NumericVector f, NumericVector f1, NumericVector f2, double N, int am_code, int sparse, int transform_code);
RcppExport SEXP wordspace_CPP_dsm_score_sparse(SEXP nrSEXP, SEXP ncSEXP, SEXP pSEXP, SEXP row_ofSEXP, SEXP fSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP NSEXP, SEXP am_codeSEXP, SEXP sparseSEXP, SEXP transform_codeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nr(nrSEXP );
        Rcpp::traits::input_parameter< int >::type nc(ncSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type row_of(row_ofSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type f1(f1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type f2(f2SEXP );
        Rcpp::traits::input_parameter< double >::type N(NSEXP );
        Rcpp::traits::input_parameter< int >::type am_code(am_codeSEXP );
        Rcpp::traits::input_parameter< int >::type sparse(sparseSEXP );
        Rcpp::traits::input_parameter< int >::type transform_code(transform_codeSEXP );
        NumericVector __result = CPP_dsm_score_sparse(nr, nc, p, row_of, f, f1, f2, N, am_code, sparse, transform_code);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
