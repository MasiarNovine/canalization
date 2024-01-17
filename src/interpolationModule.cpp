// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

// DEPRECATED, MARKED TO BE REMOVED
/**
* @brief Get the indices of age values in a reference vector.
*
* @param age The age value to find indices for.
* @param refAge The reference vector of age values.
* @return A vector with two indices, indicating the position of age in refAge.
*/
// [[Rcpp::export]]
arma::vec getIndicesCpp(float age, const arma::vec& refAge) {
    // Initialize the result vector
    arma::vec indices(2, arma::fill::zeros);

    // Calculate the difference between age and the first element of refAge
    float delta = age - refAge(0);

    // Check if delta is less than 0
    if (delta < 0) {
        // If delta is less than 0, set the result indices to 0 and 1
        indices(0) = 0;
        indices(1) = 1;
    } else {
        // Loop through the refAge vector until delta is less than 0 or index reaches the end
        arma::uword index = 0;
        while (delta >= 0 && index < refAge.n_elem) {
            // Calculate the difference between age and the current element of refAge
            delta = age - refAge(index);

            // If delta is greater than or equal to 0, increment the index
            if (delta >= 0) {
                index++;
            }
        }

        // Set the result indices to index - 1 and index
        indices(0) = index - 1;
        indices(1) = index;
    }

    // Return the result vector
    return indices;
}

// DEPRECATED, MARKED TO BE REMOVED
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

// POSSIBLE DEPRECATED, MARKED TO BE REMOVED
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
    float interpolatedValue = lower + (1 - proportion) * (upper - lower);
    //Rcpp::Rcout << "lower + (1 - proportion) * (upper - lower) = " << lower << " + " << (1 - proportion) << " * " << upper << " - " << lower << " = " << interpolatedValue << std::endl;
    return interpolatedValue;
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

 // [[Rcpp::export]]
float calculateCentileCpp(float zValue, float power, float median, float coeffVar) {
    // calculate the z-score
    float centile = median * pow(1 + power * coeffVar * zValue, 1 / power);

    // return the z-score
    return centile;
}

// POSSIBLY DEPRECATED, MARKED TO BE REMOVED
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
    arma::vec result(inputValues.n_elem, arma::fill::zeros);

    for (arma::uword i = 0; i < inputValues.n_elem; i++) {
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

//[[Rcpp::export]]
arma::vec findAgeIndicesCpp(const arma::vec& inputAges, const Rcpp::StringVector& inputSexes,
                            const arma::vec& referenceAges, const arma::uvec& referenceFemaleIndices,
                            const arma::uvec& referenceMaleIndices
) {
    arma::vec indices(inputAges.n_elem, arma::fill::zeros);

    if (inputAges.is_empty() || referenceAges.is_empty() ||
        referenceFemaleIndices.is_empty() || referenceMaleIndices.is_empty()) {
        return indices;
    }

    for (arma::uword i = 0; i < inputAges.n_elem; i++) {
        const arma::uvec& referenceIndices = (inputSexes(i) == "female") ? referenceFemaleIndices : referenceMaleIndices;

        if (!referenceIndices.empty()) {
            arma::uword index = 0;
            float delta = referenceAges(referenceIndices(0)) - inputAges(i);

            while (delta != 0 && index < referenceIndices.n_elem) {
                delta = referenceAges(referenceIndices(index)) - inputAges(i);
                if (delta != 0) index++;
            }

            indices(i) = index;
        }
    }

    return indices;
}

//[[Rcpp::export]]
arma::vec mapAgeValuesCpp(const arma::vec& inputAgeValues, const arma::vec& referenceAgeValues) {
    arma::vec indices(inputAgeValues.n_elem, arma::fill::zeros);

    if (inputAgeValues.is_empty() || referenceAgeValues.is_empty()) {
        return indices;
    }

    for (arma::uword i = 0; i < inputAgeValues.n_elem; i++) {
        arma::uword index = 0;
        float delta = referenceAgeValues(0) - inputAgeValues(i);

        while (delta != 0 && index < referenceAgeValues.n_elem) {
            delta = referenceAgeValues(index) - inputAgeValues(i);
            if (delta != 0) index++;
        }

        indices(i) = index;
    }

    // CAREFUL: THIS FUNCTION SPECIFIES THE INDEX FOR R!
    return indices + 1;
}

//[[Rcpp::export]]
arma::vec calculateZscoresCpp(const arma::vec& inputValues, const arma::vec& inputAges, const Rcpp::StringVector& inputSexes,
                              const arma::vec& referenceAges, const arma::vec& referenceL, const arma::vec& referenceM,
                              const arma::vec& referenceS, const arma::uvec& referenceFemaleIndices, const arma::uvec& referenceMaleIndices
) {
    arma::vec result(inputValues.n_elem, arma::fill::zeros);
    arma::vec indices = findAgeIndicesCpp(inputAges,inputSexes, referenceAges,
                                          referenceFemaleIndices, referenceMaleIndices);

    if (!indices.is_empty()) {
        for (arma::uword i = 0; i < inputValues.n_elem; i++) {
            result(i) = calculateZScoreCpp(inputValues(i), referenceL(indices(i)), referenceM(indices(i)), referenceS(indices(i)));
        }
    }

    return result;
}

//[[Rcpp::export]]
arma::vec calculateCentilesCpp(const arma::vec& zValues, const arma::vec& inputAges, const Rcpp::StringVector& inputSexes,
                               const arma::vec& referenceAges, const arma::vec& referenceL, const arma::vec& referenceM,
                               const arma::vec& referenceS, const arma::uvec& referenceFemaleIndices, const arma::uvec& referenceMaleIndices
) {
    arma::vec result(zValues.n_elem, arma::fill::zeros);
    arma::vec indices = findAgeIndicesCpp(inputAges, inputSexes, referenceAges,
                                          referenceFemaleIndices, referenceMaleIndices);

    if (indices.is_empty()) {
        for (arma::uword i = 0; i < zValues.n_elem; i++) {
            result(i) = calculateCentileCpp(zValues(i), referenceL(indices(i)), referenceM(indices(i)), referenceS(indices(i)));
        }
    }

    return result;
}

// This function is not used anymore.
arma::vec calculateCutoffCentilesCpp(const arma::vec& zValues, const arma::vec& referenceL,
                                     const arma::vec& referenceM, const arma::vec& referenceS
) {
    arma::vec result(zValues.n_elem, arma::fill::zeros);

    for (arma::uword i = 0; i < zValues.n_elem; i++) {
        result(i) = calculateCentileCpp(zValues(i), referenceL(i), referenceM(i), referenceS(i));
    }

    return result;
}

// [[Rcpp::export]]
arma::vec assignMedianBmiCpp(const arma::vec& inputValues, const arma::vec& inputAges,
                             const Rcpp::StringVector& inputSexes, const arma::vec& referenceAges,
                             const arma::vec& referenceM, const arma::uvec& referenceFemaleIndices,
                             const arma::uvec& referenceMaleIndices
) {
    arma::vec result(inputValues.n_elem, arma::fill::zeros);
    arma::vec indices = findAgeIndicesCpp(inputAges, inputSexes, referenceAges,
                                          referenceFemaleIndices, referenceMaleIndices);

    if (!indices.is_empty()) {
        for (arma::uword i = 0; i < inputValues.n_elem; i++) {
            result(i) = referenceM(indices(i));
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
    for (arma::uword i = 0; i < ageBounds.n_rows; ++i) {
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
    arma::uword index = 0;

    // Loop through all subject IDs
    for (arma::uword i = 0; i < subjectIds.length(); i++) {
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

/**
 * Finds the indices of elements in a vector that match the elements in a reference vector.
 *
 * @param ages the vector of ages to search through
 * @param referenceAges the reference vector of ages to match against
 *
 * @return a vector of indices indicating the positions of matching elements in the reference vector
 *
 * @throws None
 */
// [[Rcpp::export]]
arma::vec findIndicesCpp(const arma::vec& ages, const arma::vec& referenceAges) {
    arma::vec indices(ages.n_elem, arma::fill::zeros);

    if (ages.is_empty() || referenceAges.is_empty()) {
        return indices;
    }

    for (arma::uword i = 0; i < ages.n_elem; i++) {
        arma::uword index = 0;
        float delta = referenceAges(0) - ages(i);

        while (delta != 0 && index < referenceAges.n_elem) {
            delta = referenceAges(index) - ages(i);
            if (delta != 0) index++;
        }

        indices(i) = index;
    }

    // Specifies index for R
    return indices + 1;
}

// Rcpp module definition
RCPP_MODULE(interpolation_module) {
    Rcpp::function("getIndicesCpp", &getIndicesCpp);
    Rcpp::function("findIndicesCpp", &findIndicesCpp);
    Rcpp::function("getProportionCpp", &getProportionCpp);
    Rcpp::function("getAgeCohortsCpp", &getAgeCohortsCpp);
    Rcpp::function("getUniqueIdsCpp", &getUniqueIdsCpp);
    Rcpp::function("interpolateCpp", &interpolateCpp);
    Rcpp::function("findAgeIndicesCpp", &findAgeIndicesCpp);
    Rcpp::function("mapAgeValuesCpp", &mapAgeValuesCpp);
    Rcpp::function("calculateZScoreCpp", &calculateZScoreCpp);
    Rcpp::function("calculateZscoresCpp", &calculateZscoresCpp);
    Rcpp::function("calculateCentileCpp", &calculateCentileCpp);
    Rcpp::function("calculateCentilesCpp", &calculateCentilesCpp);
    Rcpp::function("interpolateToZscoreCpp", &interpolateToZscoreCpp);
    Rcpp::function("assignMedianBmiCpp", &assignMedianBmiCpp);
}