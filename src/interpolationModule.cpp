#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
/**
* @brief Get the indices of age values in a reference vector.
*
* @param age The age value to find indices for.
* @param ref_age The reference vector of age values.
* @return A vector with two indices, indicating the position of age in ref_age.
*/
// [[Rcpp::export]]
arma::vec getIndicesCpp(float age, const arma::vec& ref_age) {
    // Initialize the result vector
    arma::vec res(2, arma::fill::zeros);

    // Initialize the index variable
    int index = 0;

    // Calculate the difference between age and the first element of ref_age
    float delta = age - ref_age(0);

    // Check if delta is less than 0
    if (delta < 0) {
        // If delta is less than 0, set the result indices to 0 and 1
        res(0) = index;
        res(1) = index + 1;
    } else {
        // Loop through the ref_age vector until delta is less than 0 or index reaches the end
        while (delta >= 0 && index < ref_age.n_rows - 1) {
            // Calculate the difference between age and the current element of ref_age
            delta = age - ref_age[index];

            // If delta is greater than or equal to 0, increment the index
            if (delta >= 0) {
                index++;
            }
        }

        // Set the result indices to index - 1 and index
        res(0) = index - 1;
        res(1) = index;
    }

    // Return the result vector
    return res;
}

/**
 * @brief Calculate the proportion of age relative to a reference age range
 *
 * @param age The age value to calculate the proportion for
 * @param refAge The reference age vector
 * @param indices The indices of the reference age range
 *
 * @return The proportion of age relative to the reference age range
 */
float getProportionCpp(float age, const arma::vec& refAge, const arma::vec& indices) {
    // Calculate the difference between the upper and lower bounds of the reference age range
    float deltaRefAge = refAge(indices(1)) - refAge(indices(0));

    // Calculate the difference between the upper bound of the reference age range and the age value
    float deltaAgeToRef = refAge(indices(1)) - age;

    // Calculate and return the proportion of age relative to the reference age range
    return deltaAgeToRef / deltaRefAge;
}

 /**
  * @brief Calculate the z-score transformation of a given value.
  *
  * @param value The value to transform.
  * @param power The power to raise value to.
  * @param median The median of the distribution.
  * @param coeffVar The coefficient of variation (normalized standard deviation) of the distribution.
  *
  * @return The z-score transformed value.
  */
 // [[Rcpp::export]]
 float calculateZScoreCpp(float value, float power, float median, float coeffVar) {
     // calculate the z-score
     float zScore = (pow(value / median, power) - 1) / (power * coeffVar);

     // return the z-score
     return zScore;
 }

/**
 * @brief Interpolates a value between two given values based on a proportion.
 *
 * @param proportion The proportion between the lower and upper values.
 * @param lower The lower value.
 * @param upper The upper value.
 *
 * @return The interpolated value between the lower and upper values.
 */
// [[Rcpp::export]]
float interpolateCpp(float proportion, float lower, float upper) {
    // Calculate the interpolated value using the formula: interpolatedValue = upper - proportion * (upper - lower)
    float interpolatedValue = upper - proportion * (upper - lower);
    return interpolatedValue;
}

/**
 * @brief Interpolates input values to z-scores based on reference data.
 *
 * @param inputValues Input values to be interpolated.
 * @param inputAges Ages corresponding to the input values.
 * @param inputSexes Sexes corresponding to the input values.
 * @param referenceAges Ages corresponding to the reference data.
 * @param referenceL Reference L values for interpolation.
 * @param referenceM Reference M values for interpolation.
 * @param referenceS Reference S values for interpolation.
 * @param referenceFemaleIndices Indices of female reference data.
 * @param referenceMaleIndices Indices of male reference data.
 *
 * @return Interpolated values in z-scores.
 */
// [[Rcpp::export]]
arma::vec interpolateToZscoreCpp(const arma::vec& inputValues, const arma::vec& inputAges, const Rcpp::StringVector& inputSexes,
                                 const arma::vec& referenceAges, const arma::vec& referenceL, const arma::vec& referenceM,
                                 const arma::vec& referenceS, const arma::uvec& referenceFemaleIndices, const arma::uvec& referenceMaleIndices
) {
    arma::vec result(inputValues.n_rows, arma::fill::zeros);

    for (int i = 0; i < inputValues.n_rows; i++) {
        const arma::uvec& referenceIndices = (inputSexes(i) == "female") ? referenceFemaleIndices : referenceMaleIndices;
        if (!referenceIndices.is_empty()) {
            const arma::vec& indices = getIndicesCpp(inputAges(i), referenceAges(referenceIndices));
            if (!indices.is_empty()) {
                const float proportion = getProportionCpp(inputAges(i), referenceAges(referenceIndices), indices);
                const float interpolL = interpolateCpp(proportion, referenceL(indices(0)), referenceL(indices(1)));
                const float interpolM = interpolateCpp(proportion, referenceM(indices(0)), referenceM(indices(1)));
                const float interpolS = interpolateCpp(proportion, referenceS(indices(0)), referenceS(indices(1)));
                result(i) = calculateZScoreCpp(inputValues(i), interpolL, interpolM, interpolS);
            }
        }
    }

    return result;
}

/**
 * @brief Returns age cohorts based on age bounds.
 *
 * @param ages A vector of ages.
 * @param ageBounds A matrix of age bounds.
 *
 * @return A list containing age lower and age upper vectors.
 */
// [[Rcpp::export]]
Rcpp::List getAgeCohortsCpp(const arma::vec& ages, const arma::mat& ageBounds) {
    // Initialize ageLower and ageUpper vectors with zeros
    arma::vec ageLower(ages.n_rows, arma::fill::zeros);
    arma::vec ageUpper(ages.n_rows, arma::fill::zeros);

    // Iterate over age bounds
    for (int i = 0; i < ageBounds.n_rows; ++i) {
        // Check if age falls within the bounds
        arma::uvec indices = arma::find(ages >= ageBounds(i, 0) && ages < ageBounds(i, 1));

        // Assign age bounds to ageLower and ageUpper vectors
        ageLower.elem(indices).fill(ageBounds(i, 0));
        ageUpper.elem(indices).fill(ageBounds(i, 1));
    }

    // Return a list with age lower and age upper vectors
    return Rcpp::List::create(
        Rcpp::Named("age_lower") = ageLower,
        Rcpp::Named("age_upper") = ageUpper
    );
}

/**
 * @brief Function to get unique IDs from a vector of subject IDs.
 *
 * @param subjectIds A vector of subject IDs.
 *
 * @return A vector of unique IDs.
 */
// [[Rcpp::export]]
Rcpp::CharacterVector getUniqueIdsCpp(Rcpp::CharacterVector subjectIds) {
    // Create a temporary vector to store unique IDs
    Rcpp::CharacterVector uniqueIdsTmp(subjectIds.length() / 4);

    // Set the first element of the unique IDs vector to the first subject ID
    uniqueIdsTmp[0] = subjectIds[0];

    // Initialize the index variable
    int index = 0;

    // Loop through all subject IDs
    for (int i = 0; i < subjectIds.length(); i++) {
        // Check if the current subject ID is already in the unique IDs vector
        while (uniqueIdsTmp(index) != subjectIds(i)) {
            // Increment the index variable
            index++;
            // Add the subject ID to the unique IDs vector
            uniqueIdsTmp(index) = subjectIds(i);
        }
    }

    // Create a vector of indices from 0 to the final index
    Rcpp::IntegerVector indices = Rcpp::seq(0, index);

    // Return the subset of unique IDs using the indices vector
    return uniqueIdsTmp[indices];
}

// Rcpp module definition
RCPP_MODULE(interpolation_module) {
    Rcpp::function("getIndicesCpp", &getIndicesCpp);
    Rcpp::function("getProportionCpp", &getProportionCpp);
    Rcpp::function("getAgeCohortsCpp", &getAgeCohortsCpp);
    Rcpp::function("getUniqueIdsCpp", &getUniqueIdsCpp);
    Rcpp::function("interpolateCpp", &interpolateCpp);
    Rcpp::function("calculateZScoreCpp", &calculateZScoreCpp);
    Rcpp::function("interpolateToZscoreCpp", &interpolateToZscoreCpp);
}