// [[Rcpp::plugins("cpp20")]]

//sum.cpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double rcpp_sum(NumericVector v){
    double sum{0.0};
    for (const auto &x : v){
        sum += x;
    }
    return sum;
}

// [[Rcpp::export]]
auto rcpp_lambda_1(NumericVector v, double A) {
    NumericVector res =
        sapply(v, [&](double x){return A * x;});
    return res;
}

// [[Rcpp::export]]
DataFrame rcpp_df() {
    // Creating vector v
    NumericVector v = {1, 2};

    // Creating DataFrame df
    DataFrame df = DataFrame::create(Named("V1") = v,         // simple assign
                                     Named("V2") = clone(v)); // using clone()

    // Changing vector v
    v = v * 2;

    return df;
}

// [[Rcpp::export]]
DataFrame foo(DataFrame x) {
    return(x);
}

// [[Rcpp::export]]
DataFrame AddNewCol(const DataFrame& df, std::string new_var) {
    NumericVector vec_x = df["x"];
    NumericVector vec_y = df["y"];
    df[new_var] = vec_x * Rcpp::pow(vec_y, 2);
    return df;
}

// df <- data.frame(x = runif(10), y = runif(10))
// AddNewCol( df ,"result")

// [[Rcpp::export]]
double get_time(double shape, double scale) {
    return R::rgamma(shape, scale);
}

