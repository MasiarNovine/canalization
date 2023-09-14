// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::export]]
Rcpp::CharacterVector getUniqueIdsCpp(Rcpp::CharacterVector subjectIds) {
    // Allocate an vector with sufficient size
    // The used dataset contains longitudinal data with multiple 
    // data points for each subject, so taking a fourth of the initial size of the 
    // dataset should be save.
    Rcpp::CharacterVector uniqueIdsTmp = Rcpp::rep(Rcpp::CharacterVector("NA"), subjectIds.length() / 4);

    // Include the first element
    uniqueIdsTmp[0] = subjectIds[0];
    
    // This index tells us, what the last relevant entry will be
    int index = 0;

    // Get all the unique subject Ids
    for (int i=0; i < subjectIds.length(); i++) {
        while(uniqueIdsTmp(index) != subjectIds(i)) {
            index++;
            uniqueIdsTmp(index) = subjectIds(i);
        }
    }

    Rcpp::IntegerVector indices = Rcpp::seq(0, index);

    // Extract only relevant rows
    return uniqueIdsTmp[indices];
}