// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// binary_search
IntegerVector binary_search(IntegerVector breaks, IntegerVector search, bool left_index);
RcppExport SEXP epic_binary_search(SEXP breaksSEXP, SEXP searchSEXP, SEXP left_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type search(searchSEXP);
    Rcpp::traits::input_parameter< bool >::type left_index(left_indexSEXP);
    __result = Rcpp::wrap(binary_search(breaks, search, left_index));
    return __result;
END_RCPP
}
// extract_sites
IntegerVector extract_sites(IntegerVector start, IntegerVector end, IntegerVector site, bool return_index, int min_sites);
RcppExport SEXP epic_extract_sites(SEXP startSEXP, SEXP endSEXP, SEXP siteSEXP, SEXP return_indexSEXP, SEXP min_sitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type end(endSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type site(siteSEXP);
    Rcpp::traits::input_parameter< bool >::type return_index(return_indexSEXP);
    Rcpp::traits::input_parameter< int >::type min_sites(min_sitesSEXP);
    __result = Rcpp::wrap(extract_sites(start, end, site, return_index, min_sites));
    return __result;
END_RCPP
}
// rowWhichMax
IntegerVector rowWhichMax(NumericMatrix m);
RcppExport SEXP epic_rowWhichMax(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    __result = Rcpp::wrap(rowWhichMax(m));
    return __result;
END_RCPP
}
