// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <math.h>

// Builds the age cohorts based on given age bounds with min and max values in each row of the matrix.
// [[Rcpp::export]]
Rcpp::List getAgeCohortsCpp(const arma::vec& ages, const arma::mat& age_bounds
) {
    // Allocate a vector, which includes the assignments, to which age cohorts its belongs
    arma::vec age_lower(ages.n_rows, arma::fill::zeros);
    arma::vec age_upper(ages.n_rows, arma::fill::zeros);

    // Generate as many cohorts, as rows in the age_bounds matrix
    for (int i=0; i < age_bounds.n_rows; i++) {
        // Run through all subjects and collect the relevant ones
        for (int j=0; j < ages.n_rows; j++) {
            if (ages(j) >= age_bounds(i, 0) & ages(j) < age_bounds(i, 1)) {
                // Assign the values to the respective subject id
                age_lower(j) = age_bounds(i, 0);
                age_upper(j) = age_bounds(i, 1);
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("age_lower") = age_lower, 
        Rcpp::Named("age_upper") = age_upper
    );
}