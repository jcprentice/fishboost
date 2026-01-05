#include <Rcpp.h>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix make_grm_loop(DataFrame pedigree) {
    int N = pedigree.nrows();

    NumericMatrix GRM(N, N);
    
    // NumericVector id = pedigree["id"];
    NumericVector sire = pedigree["sire"];
    NumericVector dam = pedigree["dam"];
    
    int const nan = -2147483648;
            
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // note -1 since ids labelled 1:N but C uses 0:N-1 indices
            int sire_i = sire[i] - 1;
            int sire_j = sire[j] - 1;
            int dam_i = dam[i] - 1;
            int dam_j = dam[j] - 1;
            
            if (i == j) {
                // Same individual
                GRM(i, j) = 1.0;
            } else if (sire_i == nan && sire_j == nan) {
                // R_IsNA
                // both unrelated parents
            } else if (sire_i == nan || sire_j == nan) {
                // exactly one is a parent
                if (sire_i == j || dam_i == j || sire_j == i || dam_j == i) {
                    // parent and child are related
                    GRM(i, j) = 0.5;
                }
                // parent and child are unrelated
            } else if (sire_i == sire_j && dam_i == dam_j) {
                // Full sib
                GRM(i, j) = 0.5;
            } else if (sire_i == sire_j || dam_i == dam_j) {
                // Half sib
                GRM(i, j) = 0.25;
            }
            // else Unrelated, already 0
        }
    }
    
    return GRM;
}