// [[Rcpp::depends(RcppArmadillo)]] 
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
float z_score_cpp(float y, float l, float m, float s) {
    float z = (pow(y / m, l) - 1) / (l * s);

    return(z);
}

// [[Rcpp::export]]
arma::vec get_indices_cpp(float age, const arma::vec& ref_age) {
    arma::vec res(2, arma::fill::zeros);
    int index = 0;
    float delta = age - ref_age(0);

    while (delta >= 0) {
        delta = age - ref_age[index];
        index++;
    }

    if (age < ref_age(0)) {
        res(0) = index;
        res(1) = index + 1;
    } else {
        res(0) = index - 1;
        res(1) = index;
    };

    return res;
}

// [[Rcpp::export]]
float get_proportion_cpp(float age, const arma::vec& ref_age, const arma::vec& indices) {
    float delta_ref_age = ref_age(indices(1)) - ref_age(indices(0));
    float delta_age_to_ref = ref_age(indices(1)) - age;

    return delta_age_to_ref / delta_ref_age;
}

// [[Rcpp::export]]
float lint_cpp(float prop, float low, float up) {
    return up - (prop * (up - low));
}

// [[Rcpp::export]]
float interpolate_to_z_score_cpp(float value, float age, const arma::vec& indices, float prop,
                               const arma::vec& ref_l, const arma::vec& ref_m, const arma::vec& ref_s
) {
    float interpolL = lint_cpp(prop, ref_l(indices(0)), ref_l(indices(1)));
    float interpolM = lint_cpp(prop, ref_m(indices(0)), ref_m(indices(1)));
    float interpolS = lint_cpp(prop, ref_s(indices(0)), ref_s(indices(1)));

    return z_score_cpp(value, interpolL, interpolM, interpolS);
}

// For given raw measurements convert them all to z-scores.
// [[Rcpp::export]]
arma::vec interpolate_to_z_score_vector_cpp(const arma::vec& values, const arma::vec& ages, const arma::vec& ref_age,
                                          const arma::vec& ref_l, const arma::vec& ref_m, const arma::vec& ref_s
) {
    arma::vec res(values.n_rows, arma::fill::zeros);
    arma::vec indices;
    float prop;
    float z;

    for (int i=0; i < values.n_rows; i++) {
        indices = get_indices_cpp(ages(i), ref_age);
        prop = get_proportion_cpp(ages(i), ref_age, indices);
        res(i) = interpolate_to_z_score_cpp(values(i), ages(i), indices, prop, ref_l, ref_m, ref_s);
    }

    return res;
}

// [[Rcpp::export]]
Rcpp::List get_age_cohorts_cpp(const arma::vec& ages, const arma::mat& age_bounds
) {
    // Allocate a vector, which includes the assignments, to which age cohorts its belongs
    arma::vec age_lower(ages.n_rows, arma::fill::zeros);
    arma::vec age_upper(ages.n_rows, arma::fill::zeros);

    // Generate as many cohorts, as rows in the age_bounds matrix
    for (int i=0; i < age_bounds.n_rows; i++) {
        // Run through all observations and collect the relevant ones
        for (int j=0; j < ages.n_rows; j++) {
            if (ages(j) >= age_bounds(i, 0) & ages(j) < age_bounds(i, 1)) {
                // Assign the values to the respective observation
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

// [[Rcpp::export]]
Rcpp::CharacterVector get_unique_ids_cpp(Rcpp::CharacterVector subject_ids) {
    // Allocate an vector with sufficient size
    // The used dataset contains longitudinal data with multiple 
    // data points for each subject, so taking a fourth of the initial size of the 
    // dataset should be save.
    Rcpp::CharacterVector unique_ids_tmp = Rcpp::rep(Rcpp::CharacterVector("NA"), subject_ids.length() / 4);

    // Include the first element
    unique_ids_tmp[0] = subject_ids[0];
    
    // This index tells us, what the last relevant entry will be
    int index = 0;

    // Get all the unique subject Ids
    for (int i=0; i < subject_ids.length(); i++) {
        while(unique_ids_tmp(index) != subject_ids(i)) {
            index++;
            unique_ids_tmp(index) = subject_ids(i);
        }
    }

    Rcpp::IntegerVector indices = Rcpp::seq(0, index);

    // Extract only relevant rows
    return unique_ids_tmp[indices];
}

// Rcpp module definition
RCPP_MODULE(interpolation_module) {
    using namespace Rcpp;

    function("get_indices_cpp", &get_indices_cpp);
    function("get_proportion_cpp", &get_proportion_cpp);
    function("get_age_cohorts_cpp", &get_age_cohorts_cpp);
    function("get_unique_ids_cpp", &get_unique_ids_cpp);
    function("lint_cpp", &lint_cpp);
    function("z_score_cpp", &z_score_cpp);
    function("interpolate_to_z_score_cpp", &interpolate_to_z_score_cpp);
    function("interpolate_to_z_score_vector_cpp", &interpolate_to_z_score_vector_cpp);
}