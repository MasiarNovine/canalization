// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <math.h>

// General purpose function for converting raw measurements to z-scores.
// [[Rcpp::export]]
float transformToZscoreCpp (float y, float L, float M, float S) {
    float zScore = (pow(y / M, L) - 1) / (L * S);

    return(zScore);
}

// The function expects the right reference age, i.e. based on the right sex
// [[Rcpp::export]]
arma::vec getIndexToInterpolateCpp (float age, const arma::vec& refAge) {
    arma::vec res(2, arma::fill::zeros);
    int index = 0;
    float delta = age - refAge(0);

    while (delta >= 0) {
        delta = age - refAge[index];
        index++;
    }

    if (age < refAge(0)) {
        res(0) = index;
        res(1) = index + 1;
    } else {
        res(0) = index - 1;
        res(1) = index;
    };

    return res;
}

// [[Rcpp::export]]
float getProportionToInterpolateCpp (float age, const arma::vec& refAge, const arma::vec& indices) {
    float deltaRefAge = refAge(indices(1)) - refAge(indices(0));
    float deltaAgeToRef = refAge(indices(1)) - age;

    return deltaAgeToRef / deltaRefAge;
}

// [[Rcpp::export]]
float lintCpp(float prop, float low, float up) {
    return up - (prop * (up - low));
}

// [[Rcpp::export]]
float interpolateToZscoreCpp (float value, float age, const arma::vec& indices, float prop,
                              const arma::vec& refL, const arma::vec& refM, const arma::vec& refS
) {
    float interpolL = lintCpp(prop, refL(indices(0)), refL(indices(1)));
    float interpolM = lintCpp(prop, refM(indices(0)), refM(indices(1)));
    float interpolS = lintCpp(prop, refS(indices(0)), refS(indices(1)));

    return transformToZscoreCpp(value, interpolL, interpolM, interpolS);
}

// For given raw measurements convert them all to z-scores.
// [[Rcpp::export]]
arma::vec interpolateToZscoresCpp (const arma::vec& values, const arma::vec& ages, const arma::vec& refAge,
                                   const arma::vec& refL, const arma::vec& refM, const arma::vec& refS
) {
    arma::vec res(values.n_rows, arma::fill::zeros);
    arma::vec indices;
    float prop;
    float zScore;

    for (int i=0; i < values.n_rows; i++) {
        indices = getIndexToInterpolateCpp(ages(i), refAge);
        prop = getProportionToInterpolateCpp(ages(i), refAge, indices);
        res(i) = interpolateToZscoreCpp(values(i), ages(i), indices, prop, refL, refM, refS);
    }

    return res;
}

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