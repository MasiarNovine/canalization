// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <math.h>

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

float getProportionToInterpolateCpp (float age, const arma::vec& refAge, const arma::vec& indices) {
    float deltaRefAge = refAge(indices(1)) - refAge(indices(0));
    float deltaAgeToRef = refAge(indices(1)) - age;

    return deltaAgeToRef / deltaRefAge;
}

float lintCpp(float prop, float low, float up) {
    return up - (prop * (up - low));
}

float transformToZscoreCpp (float y, float L, float M, float S) {
    float zScore = (pow(y / M, L) - 1) / (L * S);

    return(zScore);
}

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