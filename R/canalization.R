#' \name{capwords}
#' \alias{capwords}
#' \title{Capitalize the first letter of each word}
#' \usage{capwords(s, strict)}
#' \description{Capitalize the first letter of each word}
#' \arguments{
#' \item{s}{A character vector}
#' \item{strict}{A logical value. Should only the first letter be capitalized? Default: FALSE. }
#' }
#' \details{
#' The function is taken from the documentation of the \code{base::chartr} function
#' and capitalizes the first letter of each word.
#' }
#' \value{A character vector}

capwords <- function(s, strict = FALSE) {
    cap <- function(s) {
        paste(toupper(substring(s, 1, 1)), {
            s <- substring(s, 2);
            # strict = TRUE` capitalizes the first letter, but only of the first word
            if (strict) tolower(s) else s
        }, sep = "", collapse = " " )
    }
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#' \name{capfactorDT}
#' \alias{capfactorDT}
#' \title{Capitalize the first letter of each word, but keeping factors as factors}
#' \usage{capfactorDT(dt, colnames)}
#' \description{Capitalize the first letter of each word, but keeping factors as factors}
#' \arguments{
#'   \item{dt}{The data.table containing the columns.}
#'   \item{colnames}{A character vector with the names for factor columns.}
#' }
#' \details{
#' This function uses the \code{capwords} function and capitalizes the first letter of
#' each word in the given columns of a data.table based on the given column names,
#' but will keep factors as factors.
#' }
#' \value{The initial data.table with the updated column}
#' \seealso{See also \code{\link{capwords}}.}

capfactorDT <- function(dt, colnames) {
    col <- dt[[colnames]]  # Use double brackets to extract the column as a factor
    levels(col) <- capwords(levels(col))
    dt[[colnames]] <- col  # Use double brackets to assign the updated factor back to the column

    return(dt)
}

#' \name{get_indices}
#' \alias{get_indices}
#' \title{Retrieve the pair of indices for linear interpolation.}
#' \usage{get_indices(age, ref, sex)}
#' \description{Retrieve the indices for the interpolation of the
#' LMS parameters}
#' \details{
#' The references for the raw measurements only are given for discrete
#' months values, such that for ages in-between, the L, M and S parameters
#' must be interpolated. For simplicity, the values are linearly interpolated.
#' For this reason, the indices for the lower and upper bound in the
#' proximity of the the values of interest must be found by this function.
#' }
#' \arguments{
#' \item{age}{The age of the child.}
#' \item{ref}{The given reference.}
#' \item{sex}{The sex of the child.}
#' }
#' \value{
#'  The relevant indices for interpolation.
#' }

get_indices <- function(age, ref, sex=c("male", "female")) {
  sex <- match.arg(sex)
  index <- 1

  # Get the relevant subset for the given sex
  sex_ids <- which(ref["sex"] == sex)

  for (i in 1:nrow(ref)) {
    # TODO: Fails if age >> ref[i, "age"]
    delta <- age - ref[i, "age"]
    if (delta < 0) {
      # This takes care of the case, if the reference does
      # not contain values for age 0
      # We jump in the next available time period, which is
      # not quite as accurate
      # but we for example lack LMS parameters for the KiGGs data
      if (i == 1 & ref[i, "age"] != 0) {
        index <- i + 1
        break
      } else {
        index <- i
        break
      }
    }
  }

  return(c(index - 1, index))
}

#' \name{get_proportion}
#' \alias{get_proportion}
#' \title{Calculate the proportion factor for precise interpolation.}
#' \usage{get_proportion(age, ref, indices)}
#' \description{
#'  Get the relative proportion, which is necessary for interpolation.
#' }
#' \details{
#' The references for the raw measurements only are given for discrete
#' months values, such that for ages in-between, the L, M and S parameters
#' must be interpolated. For simplicity, the values are linearly interpolated.
#' For this reason, the indices for the lower and upper bound in the proximity
#' of the the values of interest must be found by this function.
#' }
#' \arguments{
#'  \item{age}{The age of the child.}
#'  \item{ref}{The reference as a data.frame.}
#'  \item{indices}{The vector of indices.}
#' }
#' \value{
#'  The proportion of the difference between the lower / upper reference age bound
#'  and the difference of the the age of the child to the upper reference age bound.
#' }
#' \seealso{See also \code{\link{get_indices}}}

get_proportion <- function(age, ref, indices) {
  delta_ref_age <- ref[indices[2], "age"] - ref[indices[1], "age"]
  delta_age_to_ref <- ref[indices[2], "age"] - age

  return(delta_age_to_ref / delta_ref_age)
}

#' \name{interpolate}
#' \alias{interpolate}
#' \title{Linear interpolation between an upper and lower value.}
#' \description{Perform linear interpolation based on the lower / upper
#' reference age bound.}
#' \arguments{
#'  \item{prop}{Proportion of the difference value should be used.}
#'  \item{low}{Lower value based on the lower reference age bound.}
#'  \item{up}{Upper value based on the upper reference age bound.}
#' }
#' \details{
#' The basic functionality to perform linear interpolation using the
#' proportion of the upper and lower age bound,
#' relevant for the respective age of a child.
#' }
#' \value{A 'numeric' value, representing the linear interpolation.}
#' \seealso{See also \code{\link{get_indices}}, \code{\link{get_proportion}},
#' \code{\link{interpolate_to_zscore}}}

interpolate <- function(prop, low, up) {
  return(up - (prop * (up - low)))
}

#' \name{interpolate_to_zscore}
#' \alias{interpolate_to_zscore}
#' \title{Linear interpolation of LMS parameters.}
#' \description{
#'  Interpolation of raw values of either weight, height or BMI to z-scores.
#' }
#' \arguments{
#'  \item{value}{Value to interpolate.}
#'  \item{age}{The age of the child.}
#'  \item{ref}{The reference data.frame.}
#'  \item{type}{The type of the value.}
#' }
#' \details{
#' For getting z-scores for values for ages, not covered by references tables,
#' linear transformation for ages in close proximity are used.
#' }
#' \value{
#' A 'vector' of z-score transforms.
#' }
#' \seealso{See also \code{\link{get_indices}}, \code{\link{get_proportion}}, \code{\link{interpolate}}.}

interpolate_to_zscore <- function(value, age, ref, type=c("Weight", "Height", "BMI")) {
  .check_measure(type)
  ty <- tolower(type)

  indices <- get_indices(age, ref)
  prop <- get_proportion(age, ref, indices)

  interpol_l <- interpolate(prop, ref[indices[1], paste0(ty, '_l')], ref[indices[2], paste0(ty, '_l')])
  interpol_m <- interpolate(prop, ref[indices[1], paste0(ty, '_m')], ref[indices[2], paste0(ty, '_m')])
  interpol_s <- interpolate(prop, ref[indices[1], paste0(ty, '_s')], ref[indices[2], paste0(ty, '_s')])

  return(transform_to_zscore(value, interpol_l, interpol_m, interpol_s))
}

#' \name{linref}
#' \alias{linref}
#' \title{Linear interpolation of LMS parameters.}
#' \description{
#'  Interpolation of raw values of either weight, height or BMI to z-scores.
#'
#' The references for the raw measurements only are given for discrete
#' months values, such that for ages in-between, the L, M and S parameters
#' are linearly interpolated.
#'
#' For getting z-scores for values for ages, not covered by references tables,
#' linear transformation for ages in close proximity are used.
#' }
#' \arguments{
#'  \item{ref}{The reference data.frame.}
#' }

linref <- function(ref, age_range = c(0, 18), step = 0.01) {
  age <- sex <- measure <- param <- nu <- mu <- sigma <- NULL
  inter <- rbind(
    ref[, approx(x = age, y = nu, xout = seq(age_range[1], age_range[2], step)), by = list(sex, measure)][, param := "nu"],
    ref[, approx(x = age, y = mu, xout = seq(age_range[1], age_range[2], step)), by = list(sex, measure)][, param := "mu"],
    ref[, approx(x = age, y = sigma, xout = seq(age_range[1], age_range[2], step)), by = list(sex, measure)][, param := "sigma"]
  )
  inter[, param := factor(param, levels = c("nu", "mu", "sigma"))]
  setnames(inter, c("x", "y"), c("age", "value"))
  return(inter)
}

#' \name{transform_to_zscore}
#' \alias{transform_to_zscore}
#' \title{z-score transformation}
#' \usage{transform_to_zscore(y, l, m, s)}
#' \description{
#' Transform raw biometrical measurements to z-scores
#' }
#' \details{
#' Given the LMS parameters and raw measurements for weight, height or BMI,
#' the given values are transformed to z-values, based on the LMS method of
#' Cole (1991), which is commonly applied for growth curve analysis.
#' }
#' \arguments{
#'   \item{y}{Vector of values.}
#'   \item{l}{Exponent parameter of the Box-Cox function.}
#'   \item{m}{The median value.}
#'   \item{s}{The standard coefficient of variance.}
#' }
#' \value{`vector` of transformed values.}
#' \references{
#' Cole, TJ, "The LMS method for constructing normalized growth standards",
#' European journal of clinical nutrition 44, 1 (1990), pp. 45-60.
#' }

transform_to_zscore <- function(y, l, m, s) {
  return(((y / m)^l - 1) / (l * s))
}

#' \name{transform_to_centile}
#' \alias{transform_to_centile}
#' \title{Transformation of z-scores to centiles}
#' \usage{transform_to_centile(z, l, m, s)}
#' \description{
#' Transform z-scores to centiles
#' }
#' \arguments{
#'   \item{z}{Vector of z-scores.}
#'   \item{l}{Exponent parameter of the Box-Cox function.}
#'   \item{m}{The median value.}
#'   \item{s}{The standard coefficient of variance.}
#' }
#' \value{`vector` of transformed values.}
#' \references{
#' Cole, MJ, "The LMS method for constructing normalized growth standards",
#' European journal of clinical nutrition 44, 1 (1990), pp. 45-60.
#' }
#' \seealso{
#' \code{\link{transform_to_zscore}}
#' }

transform_to_centile <- function(z, l, m, s) {
  y <- m * (1 + l * s * z)^(1/l)
  return(y)
}

#' \name{CQV}
#' \alias{CQV}
#' \title{CQV}
#' \usage{CQV(x)}
#' \description{Calculates the CQV of a vector of values.}
#' \arguments{
#'   \item{x}{Vector of values.}
#' }
#' \value{`vector` of transformed values.}

CQV <- function(x) {
 Q3 <- quantile(x, probs=c(0.75));
 Q1 <- quantile(x, probs=c(0.25));
 cqv <- IQR(x) / (Q3 + Q1);
 return(100 * cqv);
}

#' \name{get_reference_names}
#' \alias{get_reference_names}
#' \title{List of available reference tables.}
#' \description{List of available reference tables.}
#' \usage{get_reference_names()}
#' \value{`vector` of available reference tables.}

get_reference_names <- function() {
  return(c("WHO", "KiGGS", "Kromeyer-Hauschild", "IOTF"))
}

#' \name{load_ref}
#' \alias{load_ref}
#' \title{Load references.}
#' \description{Loads LMS references.}
#' \usage{load_ref(ref_name = "Kromeyer-Hauschild", measure = c("Weight", "Height", "BMI", "all"))}
#' \arguments{
#'   \item{ref_name}{Name of the reference table.}
#'   \item{measure}{Vector of measures.}
#' }
#' \details{
#' Loads reference tables for measures `Weight`, `Height` and `BMI`.
#' The reference values used are based on the GAMLSS fits provided by the
#' \code{childsds} package (author: M.Vogel), such that the parameter names
#' are nu (L), mu (M) and sigma (S).
#' }
#' \value{`data.table` with nu, mu and sigma parameters.}
#' \references{
#' TJ Cole, "The LMS method for constructing normalized growth standards", European journal of clinical nutrition 44, 1 (1990), pp. 45-60,
#' Stolzenberg, H., H. Kahl, and K. E. Bergmann (May 2007). “Körpermaße bei Kindern und Jugendlichen in Deutschland”. In: Bundesgesundheitsblatt - Gesundheitsforschung - Gesundheitsschutz 50.5, pp. 659–669. issn: 1437-1588. url: https://doi.org/10.1007/s00103-007-0227-5,
#' Kromeyer-Hauschild, K. et al. (Aug. 2001). “Perzentile für den Body-mass-Index für das Kindes- und Jugendalter unter Heranziehung verschiedener deutscher Stichproben”. In: Monatsschrift Kinderheilkunde 149.8, pp. 807–818. issn: 1433-0474. url: https://doi.org/10.1007/s001120170107.
#' Rigby, R.A. and D.M. Stasinopoulos (2005). “Generalized additive models for location, scale and shape.” In: Applied Statistics 54, pp. 507–554.
#' }

load_ref <- function(ref_name = c("Kromeyer-Hauschild", "KiGGS", "WHO", "IOTF"),
                     measure = c("all", "Weight", "Height", "BMI")) {
  ref_name <- tolower(match.arg(ref_name))
  measure <- tolower(match.arg(measure))
  ref <- switch(ref_name,
    "kromeyer-hauschild" = .map_ref(ref_name, measure),
    "kiggs" = .map_ref(ref_name, measure),
    "who" = .map_ref(ref_name, measure),
    "iotf" = .map_ref(ref_name, measure)
  )
  return(ref)
}

#' \name{lincutoffs}
#' \alias{lincutoffs}
#' \usage{lincutoffs(ref_name, cutoffs, groups, years, step)}
#' \title{Interpolation of percentile BMI cutoffs}
#' \description{Interpolate percentile BMI cutoffs for a given reference.}
#' \details{
#'  The function linearly interpolates percentile cutoffs for a given reference based on the
#'  International Obesity Task Force (IOTF) cutoffs. The reference is interpolated to the age range
#'  defined by `years` and `step`. The function returns a `data.table` with the
#'  percentile cutoffs.
#' }
#' \arguments{
#'   \item{ref_name}{Reference name. Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".}
#'   \item{cutoffs}{Cutoffs.}
#'   \item{groups}{Names of the weight groups.}
#'   \item{years}{Age range in years.}
#'   \item{step}{Age step.}
#' }
#' \value{`data.table` of percentile cutoffs.}
#' \references{
#' Cole, MJ, "The LMS method for constructing normalized growth standards",
#' European journal of clinical nutrition 44, 1 (1990), pp. 45-60.
#' }

lincutoffs <- function(ref_name = get_reference_names(),
                       cutoffs = c(18.5, 25.0, 30.0),
                       groups = c("underweight", "overweight", "obese"),
                       years = c(0, 18),
                       step = 0.01
) {
  age <- sex <- bmi_l <- bmi_m <- bmi_s <- z <- bmi_centile <- status <- patterns <- NULL

  # Check, whether the reference is given
  .check_ref_name(ref_name)

  # Interpolated LMS parameter table
  ref_interpl <- linref(ref_name = ref_name, age_range = years, step = step)

  # Calculate interpolated percentile cutoffs based on cutoffs at age 18 years
  prms <- ref_interpl[age == 18, list(bmi_l, bmi_m, bmi_s), by = sex]
  cts <- prms[, list(z = transform_to_zscore(cutoffs, l = bmi_l, m = bmi_m, s = bmi_s)), by = sex]

  # Repeat the interpolated parameters as many times as cutoffs
  # for each sex
  tmp <- vector("list", length = length(ref_interpl[, levels(sex)]))
  names(tmp) <- ref_interpl[, levels(sex)]

  for (s in ref_interpl[, levels(sex)]) {
    tmp[[s]] <- ref_interpl[sex == s, list(
        age = rep(age, times = length(cutoffs)),
        sex = rep(sex, times = length(cutoffs)),
        bmi_l = rep(bmi_l, times = length(cutoffs)),
        bmi_m = rep(bmi_m, times = length(cutoffs)),
        bmi_s = rep(bmi_s, times = length(cutoffs)),
        z = cts[sex == s, rep(z, each = ref_interpl[sex == s, .N])],
        status = rep(groups, each = ref_interpl[sex == s, .N]))]
  }

  # Collect sex specific data
  res <- rbindlist(tmp)

  # Convert status to factor
  res[, status := factor(status, levels = groups)]

  # Transform to bmi centiles
  res <- res[, bmi_centile := transform_to_centile(z, bmi_l, bmi_m, bmi_s)]

  return(res)
}

#' \name{nutritional_status_bmi}
#' \alias{nutritional_status_bmi}
#' \title{Weight status assignment}
#' \description{Assigns a weight status for each observation in the data.}
#' \details{
#' Based on given BMI cutoffs, the adiposity status for each observation is determined.
#' The argument `cutoffs` must be a named list of 2-element vectors:
#' The first element denotes the lower (inclusive) bound, while the second element
#' denotes the upper (exclusive) bound.
#' }
#' \arguments{
#'   \item{cresc_data}{A data.table containing CrescNet data.}
#'   \item{ref_name}{Reference name. Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".}
#'   \item{bounds}{A named list with BMI cutoffs given as a 2-element vector. Default: list(underweight = c(NA, 18.5), normal = c(18.5, 25), overweight = c(25, 30), obese = c(30, NA))}
#'   \item{years}{Age range in years. Default: c(0, 18)}
#'   \item{step}{Age step. Default: 0.01}
#' }
#' \value{A data.table containing an extra column with the assigned weight status for each observation.}

# TODO: Put this function and the other one for BMI-SDS together
nutritional_status_bmi <- function(cresc_data,
                               ref_name = get_reference_names(),
                               bounds = list(underweight = c(NA, 18.5), normal = c(18.5, 25), overweight = c(25, 30), obese = c(30, NA)),
                               years = c(0, 18),
                               step = 0.01
) {
  bmi <- sex <- tmp_id <- age <- status <- bmi_centile <- idxlow <- idxup <- cutlow <- cutup <- st <- s <- NULL

  # Collapse the cutoffs to a vector
  cutoffs <- unlist(bounds)

  # Interpolation reference cutoffs
  cuttab <- lincutoffs(ref_name = ref_name,
                                cutoffs = cutoffs,
                                groups = names(cutoffs),
                                years = years,
                                step = step)[, list(age, sex, status, bmi_centile)]

  # Temporary id column
  cuttab[, tmp_id := seq(1, .N, 1)]

  # Alocation of the status
  cresc_data[, status := ""]
  cresc_data[, tmp_id := seq(1, .N, 1)]
  setkey(cresc_data, age)       # Sort DT in place
  cresc_data[, .SD, by = list(age)]

  # Split by sexes
  tmp <- list(female = cresc_data[sex == levels(sex)[1], ],
              male = cresc_data[sex == levels(sex)[2], ])

  # Go trough all adiposity status
  for (st in names(bounds)) {
    for (s in cresc_data[, levels(sex)]) {

      # Get the lower and upper bound
      cutlow <- cuttab[sex == s & status == paste0(st, "1"), ]
      cutup <- cuttab[sex == s & status == paste0(st, "2"), ]

      # Index for relevant ages and sex
      idxlow <- cresc_data[sex == s, findInterval(age, cutlow[, age])]
      idxup <- cresc_data[sex == s, findInterval(age, cutup[, age])]

      # Get the BMI
      bmilow <- if (cutlow[idxlow, !all(is.na(bmi_centile))]) cutlow[idxlow, bmi_centile] else 0
      bmiup <- if (cutup[idxup, !all(is.na(bmi_centile))]) cutup[idxup, bmi_centile] else .Machine$integer.max

      # Update data.table
      tmp[[s]][bmi >= bmilow & bmi < bmiup, status := st]
    }
  }

  # Tidy up
  res <- rbind(tmp[[cresc_data[, levels(sex)][1]]], tmp[[cresc_data[, levels(sex)][2]]])
  res <- res[, .SD, keyby = tmp_id]

  return(res[, tmp_id := NULL])
}

#' \name{nutritional_status_bmi_sds}
#' \alias{nutritional_status_bmi_sds}
#' \usage{nutritional_status_bmi_sds(cresc_data, cutoffs)}
#' \title{Weight status assignment}
#' \description{Assigns a weight status for each observation in the data.}
#' \details{
#' Based on given percentile cutoffs and categories, each observation is categorized.
#' The argument `cutoffs` must be a named list of lists:
#' The first list denotes the category, while the second list includes following entries
#' with the following named entries: `lower`, `upper`, `inclower`, `incupper`, i.e. the
#' lower / upper bound (as a numeric value), whether the lower / upper bounds should be
#' inclusive or not (as a boolean), respectively. `lower` and `upper` are mandatory, while
#' the two other entries are not. If not given, the default assumes inclusive bounds.
#' }
#' \arguments{
#'   \item{cresc_data}{A data.table containing CrescNet data.}
#'   \item{cutoffs}{BMI cutoffs used.}
#' }
#' \value{A data.table containing an extra column with the assigned weight status for each observation.}

nutritional_status_bmi_sds <- function(
  cresc_data,
  cutoffs = list(
    underweight = list(lower=NA, upper=qnorm(0.05), inclower=TRUE, incupper=FALSE),
    normal = list(lower=qnorm(0.05), upper=qnorm(0.85), inclower=TRUE, incupper=FALSE),
    overweight = list(lower=qnorm(0.85), upper=qnorm(0.95), inclower=TRUE, incupper=FALSE),
    obese = list(lower=qnorm(0.95), upper=NA, inclower=TRUE, incupper=FALSE))
) {
    if (class(cutoffs) != "list") {
      stop("Argument 'cutoffs' must be a list.")
    } else if (length(cutoffs) == 0) {
      stop("Argument 'cutoffs' must contain at least one element.")
    } else if (is.null(names(cutoffs))) {
        stop("Argument 'cutoffs' must be a named list (category) with each element being a named list (list arguments) itself.")
    }

    status <- bmi_sds <- NULL
    cresc_data[, status := ""]
    nms <- names(cutoffs)

    # Find & assign end status
    # TODO: Check the bounds, e.g. the next higher status should
    # be inclusive of the lower bound and exclusive the upper bound
    # at the moment it is inclusive lower and upper bound
    for (nm in nms) {
      catg <- names(cutoffs[[nm]])

      if (!all(c("lower", "upper") %in% catg)) {
        stop("Each named list element of 'cutoffs' must be itself a named list with at least the elements 'upper' and 'lower'. See details in ?catergorize_by_endpoint.")
      }

      lower <- if (is.na(cutoffs[[nm]]$lower)) -.Machine$integer.max else cutoffs[[nm]]$lower
      upper <- if (is.na(cutoffs[[nm]]$upper)) .Machine$integer.max else cutoffs[[nm]]$upper

      if (all(c("inclower", "incupper") %in% catg)) {
        inclower <- cutoffs[[nm]]$inclower
        incupper <- cutoffs[[nm]]$incupper

        if (inclower & incupper) {
          cresc_data[between(bmi_sds, lower, upper, incbounds = TRUE), status := nm]
        } else if (inclower & !incupper) {
            cresc_data[`>=`(bmi_sds, lower) & `<`(bmi_sds, upper), status := nm]
        } else if (!inclower & incupper) {
            cresc_data[`>`(bmi_sds, lower) & `<=`(bmi_sds, upper), status := nm]
        } else cresc_data[between(bmi_sds, lower, upper, incbounds = FALSE), status := nm]
      } else cresc_data[between(bmi_sds, lower, upper, incbounds = TRUE), status := nm]
    }

    return(cresc_data[, status := factor(status, levels = nms)])
}

#' \name{assign_median_bmi}
#' \alias{assign_median_bmi}
#' \usage{assign_median_bmi(cresc_data, ref_name)}
#' \title{Assign median BMI}
#' \description{Assign median BMI for a given reference.}
#' \details{
#'  Assigning the interpolated BMI median for a given reference for
#'  age and sex adjustment.
#' }
#' \arguments{
#'   \item{cresc_data}{Input data.frame.}
#'   \item{ref_name}{Reference name. Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".}
#' }
#' \value{`data.table` with BMI median.}

assign_median_bmi <- function(cresc_data, ref_name = "Kromeyer-Hauschild") {
  # Assign global variables
  sex <- age <- bmi <- bmi_m <- NULL

  .check_ref_name(ref_name)

  # Interpolate reference
  ref_interpl <- linref(ref_name, measure = "BMI")

  # Allocate BMI to be interpolated
  cresc_data[, bmi_m := -1000];

  # Map the age the interpolated LMS parameters
  fidx <- ref_interpl[, which(sex == levels(sex)[1])] - 1
  midx <- ref_interpl[, which(sex == levels(sex)[2])] - 1

  bmimed <- assignMedianBmiCpp(
    cresc_data[, bmi], cresc_data[, age],
    cresc_data[, sex], ref_interpl[, age],
    ref_interpl[, bmi_m], fidx,
    midx
  )

  cresc_data[, bmi_m := bmimed]

  return(cresc_data)
}

#' \name{checkdup}
#' \alias{checkdup}
#' \title{Check for duplicates}
#' \usage{checkdup(cresc)}
#' \description{Check for duplicates and keep them, if there
#' are harmless otherwise convert to NA to be imputed. 'Harmless' means
#' that the values are within 1% of the mean of the duplicates.
#' }
#' \arguments{
#'   \item{cresc}{CrescNet `data.table`.}
#' }
#' \value{`data.table` with duplicates removed.}

checkdup <- function(cresc) {
  age <- subject_id  <- weight <- height <- bmi <- index <- age_group <- NULL

  # Percentage deviation
  propdev <- function(x) {
    mndff <- mean(abs(diff(x)), na.rm = TRUE)
    mn <- mean(x, na.rm = TRUE)
    return(mndff / mn)
  }
  majvote <- function(x) {
      freqx <- table(x)
      mxi <- names(freqx)[which.max(freqx)]
      return(mxi)
  }
  procdups <- function(dups, val) {
    setna <- function(x) return(NA)
    devval <- paste0("dev_", val)
    dupdev <- dups[, propdev(get(val)), by = list(age, subject_id)]
    setnames(dupdev, "V1", devval)
    dupdev[, paste0(devval) := ifelse(get(devval) > 0.01, NA, get(devval))]
    keepkeys <- dupdev[!is.na(get(devval)), list(subject_id, age)]
    nakeys <- dupdev[is.na(get(devval)), list(subject_id, age)]
    keep <- dups[subject_id %in% keepkeys[, subject_id] & age %in% keepkeys[, age],
                 mean(get(val), na.rm = TRUE),
                 by = list(subject_id, age)]
    na <- dups[subject_id %in% nakeys[, subject_id] & age %in% nakeys[, age],
               setna(get(val)),
               by = list(subject_id, age)]
    setnames(keep, "V1", val)
    setnames(na, "V1", val)
    return(rbind(keep, na))
  }
  cresc[, index := .I]
  is_dupl <- duplicated(cresc, by = c("subject_id", "age"))
  dupidx <- sort(c(cresc[is_dupl, ][, index] - 1, cresc[is_dupl, ][, index]))
  # Indices for triplies etc. identical records must be removed
  dupidx <- unique(dupidx)
  dups <- cresc[dupidx, ]
  cresc <- cresc[!dupidx, ]
  dupproc <- merge(procdups(dups, "weight"), procdups(dups, "height"), by = c("subject_id", "age"), all = TRUE)
  dupproc[, bmi := weight / (height / 100)^2]
  bascol <- c("sex", "gest_age", "type", "age_group")
  bas <- dups[subject_id %in% dupproc[, subject_id],
              unique(.SD),
              by = list(subject_id, age),
              .SDcols = bascol]
  dupproc <- merge(bas, dupproc)
  # TODO: Age_group might not always be "interval" (?)
  dupproc[, age_group := ifelse(is.na(weight) | is.na(height), "removed", "interval")]
  cresc <- rbind(cresc[, .SD, .SDcols = colnames(dupproc)], dupproc)
  cresc <- cresc[, .SD, keyby = c("subject_id", "age")]
  return(cresc)
}

#' \name{znorm}
#' \alias{znorm}
#' \title{Z-normalization}
#' \usage{znorm(cresc, ref)}
#' \description{
#' Z-transformation of each record
#' We need to find the sex-specific indices for each individual record based on the
#' interpolated reference table, which for convenience is interpolated in 0.01 years steps.
#' }
#' \arguments{
#'   \item{cresc}{CrescNet data.table.}
#'   \item{ref}{Reference data.table given in wide format.}
#' }
#' \value{`data.table` with z-normalized values.}

znorm <- function(cresc, ref) {
  age <- sex <- weight_sds <- height_sds <- bmi_sds <- measure <- mu <- nu <- sigma <- NULL
  # z-transformation of each record
  cresc[, weight_sds := rep(-1000, .N)]
  cresc[, height_sds := rep(-1000, .N)]
  cresc[, bmi_sds := rep(-1000, .N)]
  # We need to find the sex-specific indices for each individual record based on the
  # interpolated reference table, which for convenience is interpolated in 0.01 years steps.
  for (ms in c("weight", "height", "bmi")) {
    for (sx in c("female", "male")) {
      intertmp <- ref[sex == sx & measure == ms, ]
      indices <- findIndicesCpp(cresc[sex == sx, round(age, 3)], intertmp[, round(age, 3)])
      transtmp <- intertmp[indices, transform_to_zscore(cresc[sex == sx, get(ms)], l=nu, m=mu, s=sigma)]
      coln <- paste0(ms, "_sds")
      cresc[sex == sx, (coln) := transtmp]
    }
  }
  return(cresc)
}

#' \name{cutage}
#' \alias{cutage}
#' \title{Age groups}
#' \usage{cutage(x, breaks = c(-Inf, 2, 6, 11, 14, Inf),
#'               labels = c("0-2yr", "2-6yr", "6-11yr", "11-14yr", "+14yr"),
#'               right = TRUE)}
#' \description{Age groups}
#' \details{Age groups}
#' \arguments{
#'   \item{x}{Age.}
#'   \item{breaks}{Age boundaries.}
#'   \item{labels}{Age labels.}
#'   \item{right}{Whether the intervals are right-closed (default), left-open or both.}
#'  }
#' \value{
#'   `factor` with age groups.
#'  }
#' \seealso{\code{\link{cut}}}

cutage <- function(x, breaks = c(-Inf, 2, 6, 11, 14, Inf),
                   labels = c("0-2yr", "2-6yr", "6-11yr", "11-14yr", "+14yr"),
                   right = TRUE
) {
  cut(x = x, breaks = breaks, labels = labels, right = TRUE)
}

#' \name{prepare_dataset}
#' \alias{prepare_dataset}
#' \usage{prepare_dataset(file_adiposity, file_control, n_train,
#' center_age, use_days, encode_sex, short_ids, first_extra, check_dups,
#' seed, age_range)}
#' \title{Loading and preprocessing of the CrescNet dataset.}
#' \description{Prepare the CrescNet dataset.}
#' \details{The functions expects file paths to specific CrescNet data with the variables
#' 'subject_id', 'sex', 'gestational_age', 'age', 'height', 'weight', 'height_sds',
#' 'bmi', 'bmi_sds' and 'age_group'. `check_dups` will look for duplicated values,
#' if their values are indistinguishable from eachother they averaged and kept.
#' If the percentage difference is within 1% of the mean them, there are also kept.
#' Otherwise their values are set to `NA`.
#' }
#' \arguments{
#'   \item{file_adiposity}{CrescNet file path for obese labeled subjects.}
#'   \item{file_control}{CrescNet file path for control labeled subjects.}
#'   \item{n_train}{Number of samples used for model building.}
#'   \item{is_znorm}{Include z-normalization? Default: TRUE}
#'   \item{center_age}{Should age be centered? Default: FALSE}
#'   \item{use_days}{Age given in days, instead of years? Default: FALSE}
#'   \item{encode_sex}{Sex encoded as 1 (male) and -1 (female)? Default: FALSE}
#'   \item{short_ids}{Shorter subject IDs? Default: TRUE}
#'   \item{first_extra}{First measurements as separate variables? Default: FALSE}
#'   \item{check_dups}{Check for duplicates? Default: FALSE}
#'   \item{seed}{Seed used for random sampling of subjects. Default: 1}
#'   \item{age_range}{Range of ages to keep. Default: NULL}
#' }
#' \value{
#'  A data.frame contraining all variables of the CrescNet dataset.
#' }

# TODO: Include attributes for units & labels see Hmisc::label and Hmisc::units
prepare_dataset <- function(file_adiposity,
                            file_control,
                            n_train = 100,
                            is_znorm = TRUE,
                            center_age = FALSE,
                            use_days = FALSE,
                            encode_sex = FALSE,
                            short_ids = TRUE,
                            first_extra = FALSE,
                            check_dups = FALSE,
                            seed = 1,
                            age_range = NULL
) {
  subject_id <- sex <- age <- record_id <- N <- age_group <- NULL

  # Loading CrescNet data & pool data sets
  cresc_data <- rbind(data.table::fread(file_adiposity,
                                        header = TRUE,
                                        stringsAsFactors = TRUE)[, "type" := "adiposity"],
                      data.table::fread(file_control,
                                        header=TRUE,
                                        stringsAsFactors = TRUE)[, "type" := "control"])

  setnames(cresc_data, "gestational_age", "gest_age")
  setkey(cresc_data, "subject_id")
  cresc_data[, ("type") := factor(get("type"), levels = c("control", "adiposity"))]

  # Should shorter ids been used?
  if (short_ids) {
    short_ids <- sort(cresc_data[, get("subject_id")])
    levels(short_ids) <- mcg.autonames("ID_", length(levels(short_ids)))
    cresc_data[, ("subject_id") := short_ids]
  }

  # Sould sex be encoded as a continous variable over [-1, 1] from male to female?
  if (encode_sex) {
    cresc_data[, sex := unlist(lapply(sex, function(x) return(if (x == "male") 1 else -1)))]
  }

  # Should age be set in days?
  if (use_days) {
    cresc_data <- cresc_data[, age := floor(age * 365.25)]
  }

  # Should age be centered?
  if (center_age) {
    cresc_data[, age := age - (max(age) - ((max(age) - min(age)) / 2))]
  }

  # Check duplicated observations
  if (check_dups) {
    cresc_data <- checkdup(cresc_data)
  }

  if (is_znorm) {
    cresc_data[, height_sds := NULL]
    cresc_data[, bmi_sds := NULL]
    ref <- load_ref(ref_name = "Kromeyer-Hauschild", measure = "all")
    interpl <- dcast(linref(ref), measure + sex + age ~ param, value.var = "value")
    cresc_data <- znorm(cresc_data, interpl)
  }

  # Assign the number of measurements per subject
  cresc_data[, N := rep(.N, .N), by = "subject_id"]

  # Include first response measurement as baseline covariate?
  if (first_extra) {
    measure <- c("height", "weight", "bmi")
    # Get response at baseline, i.e. first measurement at approx age 0.5 years
    start <- cresc_data[get("age_group") == "start", .SD, .SDcols = c("subject_id", "age", tolower(measure))]
    setnames(start, c("age", tolower(measure)), c("age0", paste0(tolower(measure), "0")))
    setkey(start, subject_id)
    cresc_data <- merge(cresc_data[age_group != "start", ], start, by = "subject_id")
  }

  # Enumerate measurements with ID
  cresc_data[, record_id := seq_len(.N), by = subject_id]

  # Include age cohorts
  if (!is.null(age_range)) {
    cresc_data <- create_age_cohorts(cresc_data, age_range)
  }

  # Draw randomly subjects for data splitting
  set.seed(seed)
  rndidx <- sample(seq(1, length(cresc_data[, unique(subject_id)]), 1), size = n_train)
  ids <- cresc_data[, unique(subject_id)][rndidx]

  # Split data
  return(list(train = cresc_data[subject_id %in% ids, ],
              test = cresc_data[!(subject_id %in% ids), ]))
}

# DEPRECATED use cutage() instead
#' \name{create_age_cohorts}
#' \alias{create_age_cohorts}
#' \title{Create age cohorts}
#' \description{
#' Creating age cohorts.
#' }
#' \details{
#'   Based on the CrescNet data, different age cohorts are created
#'   using given boundaries. The boundaries are given in form of a `matrix`.
#' }
#' \arguments{
#'   \item{cresc_data}{A `data.frame` with CrescNet data.}
#'   \item{age_range}{A 2*n `matrix` with age boundaries.}
#' }
#' \value{
#'   `list` with as much entries as age cohorts.
#' }

create_age_cohorts <- function(cresc_data, age_range) {
  age <- NULL

  # Assign age ranges to observations
  age_cohorts <- getAgeCohortsCpp(cresc_data[, age], age_range)

  # Create levels
  y <- gl(
    nrow(age_range), 1, nrow(age_range),
    labels=paste0(age_range[, 1], "-", age_range[, 2], " years")
  )

  # Exclude those subject, which have not been assigned to a age group
  idx <- which(age_cohorts$age_lower == 0 & age_cohorts$age_upper == 0)

  # Exclude non assigned subjects
  res <- if (length(idx) != 0) {
    # Create factors for each group but exclude to assigned subjects
    age_range <- factor(paste0(age_cohorts$age_lower[-idx], "-", age_cohorts$age_upper[-idx], " years"), levels=levels(y))
    cbind(cresc_data[-idx, ], age_range=age_range)
  } else {
    # Create factors for each group
    age_range <- factor(paste0(age_cohorts$age_lower, "-", age_cohorts$age_upper, " years"), levels=levels(y))
    cbind(cresc_data, age_range=age_range)
  }

  return(res)
}

#' \name{do_simulation}
#' \alias{do_simulation}
#' \title{Run a sampling simulation to emulate CrescNet data filtering}
#' \description{
#' Running a simulation based on the general properties of the CrescNet in use.
#' }
#' \details{
#' Simulation function to investigate the effect of the stratified sampling, i.e.
#' the use of split criteria on the association of the BMI mean and standard deviation (SD).
#'
#' }
#' \arguments{
#' \item{cresc_data:}{A `data.frame` with CrescNet data.}
#' \item{n_subjects:}{Number of simulated subjects. Default: 1000.}
#' \item{k:}{Number of observations per subject. Default: 20.}
#' \item{n_pop:}{Number of observations in the population. Default: 1e5.}
#' \item{split:}{Split criterion. Default: 0.97.}
#' }
#' \value{
#' A `list` with with two data.tables: `sample` and `aggregated`.
#' }

do_simulation <- function(cresc_data, n_subjects = 1000, k = 20, n_pop = 1e5, split = 0.97) {
    subject_id <- age_group <- record_id <- bmi <- bmi_sds <- type <- mean_var <- sd_var <- mu_bmi <- var_bmi <- age_group <- N <- NULL

    # Get the BMI-SDS variability for each groups
    varsubj <- cresc_data[, list(sd=sd(bmi_sds)), by = list(subject_id, type)]

    # Collect the mean and standard deviation of the variabiliy
    varpars <- varsubj[, list(mean_var = mean(sd), sd_var = sd(sd)), by = list(type)]

    # Define a normal distribution representing the overall population with the empirical variability
    # of each subpopulation
    pop <- varpars[, list(mu_bmi = rnorm(n_pop, 0, 1), var_bmi = rnorm(n_pop, mean_var, sd_var)), by = type]
    pop[, subject_id := mcgraph::mcg.autonames("ID_", n = .N)]

    # Sample artificial subjects with k observations with mu and sigma from each distribution
    # This is too slow right now
    bmi_pop <- pop[, list(bmi_sds = rnorm(k, mu_bmi, abs(var_bmi))), by = subject_id]

    # Select based on the 97th centile individually for each study group
    pop_ob_subj <- bmi_pop[which(bmi_sds > qnorm(split, 0, 1)), unique(subject_id)]
    bmi_pop[subject_id %in% pop_ob_subj, type := "Adiposity (simulated)"]
    bmi_pop[!(subject_id %in% pop_ob_subj), type := "Control (simulated)"]

    # Sampling
    subs_ob_subj <- sample(pop_ob_subj, n_subjects)
    subs_ctrl_subj <- sample(bmi_pop[!(subject_id %in% pop_ob_subj), subject_id], n_subjects)
    sampl <- bmi_pop[subject_id %in% c(subs_ob_subj, subs_ctrl_subj), .SD]

    # Include measurement IDs
    sampl[, record_id := seq(1, .N, 1), by = subject_id]
    sampl[, age_group := ""]
    sampl[, N := .N, by = subject_id]

    # Included age groups
    sampl[record_id == 1, age_group := "start", by = subject_id]
    sampl[record_id == N, age_group := "end", by = subject_id]
    sampl[record_id > 1 & record_id < N, age_group := "interval"]

    # Categorize adiposity groups based on BMI cutoffs
    #sampl[age_group == "end",
    #      cut(bmi, breaks = c(-Inf, pnorm(), 24.9, Inf), labels = c("Underweight", "Normal", "Obese"))]

    # Summary statistics
    agg <- sampl[, list(mean = mean(bmi_sds), sd = sd(bmi_sds)), by = list(subject_id, type)]

    return(list(sample = sampl, aggregated = agg))
}

#' \name{check_model}
#' \alias{check_model}
#' \title{Checking mixed-effects model}
#' \description{Model check applying common checks for model assessment.}
#' \usage{check_model(model, cresc_data, seed)}
#' \arguments{
#'   \item{model}{A `gls` model.}
#'   \item{cresc_data}{CrescNet data.table}
#'   \item{seed}{The seed for extracting the 20 randomly drawn subjects. Default: 1}
#' }
#' \details{
#' The usually model checks are applied, i.e. (1) constant variance residuals check, (2) random distributed residuals,
#' (3) goodness of fit by plotting fitted values against the outcome variable, (4) the variogram for observing the intra-subject correlation
#' and (5) plotting fitted values against the standardized residuals.
#' }
#' \value{Series of `ggplot` plots.}

check_model <- function(model, cresc_data, seed = 1) {
  subject_id <- age <- bmi <- NULL

  .allocate(c("bmi", "fitted", "resid", "subject_id", "rndidx", "age", "weight"))
  tmp <- copy(cresc_data)

  # Define plotting settings
  blue <- "#5499C7"
  # t <- theme_bw(base_size = 12, base_family = "Fira Sans") %+replace%
  # theme(plot.title = element_text(hjust = 0.5))
  t <- theme_elegant()

  # Update initial data table
  if (class(model) == "brmsfit") {
    # The resid is simply the difference between the raw value and the fitted estimate
    fi <- fitted(model)
    tmp[, resid := model$data[, 1] - fi[, 1]]
    tmp[, fitted := fi[, 1]]
  } else {
    tmp[, resid := resid(model)]
    tmp[, fitted := fitted(model)]
  }

  set.seed(seed)
  min_size <- if (length(tmp[, unique(subject_id)]) < 20) length(tmp[, unique(subject_id)]) else 20
  rndidx <- sample(tmp[, unique(subject_id)], size = min_size)

  # ASSUMPTION 1: CONSTANT VARIANCE
  # Check, whether the expection is 0 and the variability is constant
  p1 <- ggplot(tmp[subject_id %in% rndidx, ], aes(subject_id, resid)) +
    geom_boxplot() +
    coord_flip() +
    geom_abline(slope = 0, intercept = 0, colour = "red") +
    labs(title = expression(bold("Residuals per subject")), x = "Subject id", y = "Standardized residuals") +
    t

  # ASSUMPTION 2: RANDOM DISTRIBUTION OF RESIDUALS AROUND MEAN OF 0 AND STANDARD DEVIATION OF 1
  # Check the residuals by subject
  p2 <- ggplot(tmp[subject_id %in% rndidx, ], aes(x = fitted, y = resid, group = subject_id)) +
    geom_point(colour = blue) +
    geom_line(colour = blue) +
    facet_wrap(~ subject_id) +
    labs(title = expression(bold("Deviation from fit")), x = "Fitted", y = "Standardized residuals") +
    t

  # DIAGNOSTICS PLOT: GOODNESS OF FIT FOR RAW RESPONSE VALUES
  # Fitted line against raw values
  p3 <- ggplot(tmp[subject_id %in% rndidx, ], aes(x = age, bmi)) +
    geom_point() +
    geom_line(mapping = aes(x = age, y = fitted), colour = blue) +
    facet_wrap(~ subject_id) +
    labs(title = expression(bold("Fitted vs response")), x = "Age", y = "Response") +
    t

  if (class(model) != "brmsfit") {
    # DIAGNOSTICS PLOT: NO TREND IN WITHIN-SUBJECT CORRELATIONS
    # Variogram for checking the correlations
    p4 <- ggvario(Variogram(model, resType = "n"))
  }

  # DIAGNOSTICS PLOT: FITTED VS RESIDUALS
  p5 <- ggplot(tmp, aes(x = fitted, y = resid)) +
    geom_point(colour = blue) +
    labs(title = expression(bold("Fitted vs response")), x = "Fitted", y = "Standardized residuals") +
    t

  # Print plots sequentially
  pr <- "Enter key for next plot:"
  readline(pr); print(p1)
  readline(pr); print(p2)
  readline(pr); print(p3)
  if (class(model) != "brmsfit") { readline(pr); print(p4) }
  readline(pr); print(p5)
}

#' \name{make_predictions}
#' \alias{make_predictions}
#' \title{Generate predictions}
#' \description{Uses the given GLS model to make predictions for 20 subjects.}
#' \usage{make_predictions(model, cresc_data, seed, log)}
#' \arguments{
#'   \item{model}{A `gls` model.}
#'   \item{cresc_data}{CrescNet data.table}
#'   \item{seed}{Seed for extracting the 20 randomly drawn subjects. Default: 1}
#'   \item{log}{Should the logarithm for the response be used? Default: TRUE}
#' }
#' \details{
#' The function can be used to observe, how the model fits body weights of
#' 20 randomly selected subjects of a test set. Next to a `ggplot` plot
#' with fits for each subject, the estimated weight values are appended
#' to the initilal data.table object.
#' }
#' \value{A `ggplot` object and a `data.table` with prediction estimates.}

make_predictions <- function(model, cresc_data, seed = 1, log = TRUE) {
  subject_id <- age <- predicted <- weight <- NULL

  blue <- "#5499C7"
  set.seed(seed)

  cresc_data_test <- cresc_data$test[
    subject_id %in% sample(cresc_data$test[, unique(subject_id)], size = 20),
  ]

  cresc_data_test[, predicted := predict(model, cresc_data_test)]

  g <- if (log) {
    ggplot(cresc_data_test, aes(x = age, log(weight)))
  } else {
    ggplot(cresc_data_test, aes(x = age, weight))
  }

  g <- g +
    geom_point() +
    geom_line(mapping = aes(x = age, y = predicted), colour = blue) +
    labs(title = expression(bold("Using the GLS model to predict test data"))) +
    facet_wrap(~ subject_id) +
    theme_bw()

  print(g)

  return(cresc_data_test)
}

#' \name{ggvario}
#' \alias{ggvario}
#' \title{Plot a variogram}
#' \description{Uses ggplot to plot a variogram}
#' \usage{ggvario(v)}
#' \arguments{
#'   \item{v}{A data.frame object with columns `dist` and `variog`.}
#' }
#' \details{
#' The column `dist` contains the distances, while `variog` is the the result
#' of the variogram function.
#' }
#' \value{A `ggplot` object.}

ggvario <- function(v) {
  variog <- type <- NULL

  g <- ggplot(v, aes(x = dist, y = variog))
  g <- g +
    geom_point() +
    expand_limits(y = 0) +
    geom_smooth(method = "loess", span = 1) +
    labs(title = expression(bold("Variogram")), x = "lag (years)", y = expression(gamma(lag))) +
    theme_bw(base_family = "Fira Sans")
  g <- if ("type" %in% colnames(v)) g + facet_wrap(vars(type)) else g
  return(g)
}

#' \name{theme_elegant}
#' \alias{theme_elegant}
#' \title{Custom ggplot theme}
#' \description{A customized ggplot theme}
#' \usage{
#' theme_elegant(base_size, base_family, axis.line, axis.line.x, axis.line.y,
#' legend.title, legend.position, legend.direction, legend.justification,
#' legend.box, panel.background, panel.border, panel.grid, plot.background, strip.background,
#' strip.text, strip.text.y, ...)
#' }
#' \arguments{
#'   \item{base_size}{The base font size.}
#'   \item{base_family}{The base font family.}
#'   \item{axis.line}{The line of the axis.}
#'   \item{axis.line.x}{The line of the axis in x direction.}
#'   \item{axis.line.y}{The line of the axis in y direction.}
#'   \item{legend.title}{The title of the legend.}
#'   \item{legend.position}{The position of the legend.}
#'   \item{legend.direction}{The direction of the legend.}
#'   \item{legend.justification}{The justification of the legend.}
#'   \item{legend.box}{The box of the legend.}
#'   \item{panel.background}{The background of the panel.}
#'   \item{panel.border}{The border of the panel.}
#'   \item{panel.grid}{The grid of the panel.}
#'   \item{plot.background}{The background of the plot.}
#'   \item{strip.background}{The background of the strip.}
#'   \item{strip.text}{The text of the strip.}
#'   \item{strip.text.y}{The text of the strip in y direction.}
#'   \item{...}{Other arguments passed to `theme`.}
#' }
#' \details{
#' A customized ggplot theme.
#' }
#' \value{A `ggplot` theme.}

theme_elegant <- function(base_size = 11,
                          # base_family ="Fira Sans",
                          base_family = "Helvetica",
                          #axis.line = element_blank(),#element_line(colour = "black", linewidth = 0.25),
                          axis.line = element_line(colour = "black", linewidth = 0.25),
                          axis.line.x = NULL,#element_line(colour = "black", linewidth = 0.25),
                          axis.line.y = NULL,#element_line(colour = "black", linewidth = 0.25),
                          legend.title = element_blank(),
                          legend.position = "bottom",
                          legend.direction = "horizontal",
                          legend.justification = "center",
                          legend.box = "horizontal",
                          panel.background = element_blank(),#element_rect(fill = NA, colour = NA),
                          #panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
                          panel.border = element_blank(),
                          #panel.grid = element_line(colour = "grey92", linewidth = 0.25, linetype = "dotted"),
                          panel.grid = element_blank(),
                          #plot.background = element_blank(),
                          plot.background = element_rect(fill = "white", colour = "white"),
                          strip.background = element_rect(fill = "white", colour = "white"),
                          strip.text = element_text(face ="plain",
                                                    size = rel(0.8),
                                                    hjust = 0.5,
                                                    margin = margin(0.4, 0.4, 0.4, 0.4, unit = "cm")),
                          strip.text.y = element_text(angle = 0,
                                                      margin = margin(0.4, 0.4, 0.4, 0.4, unit = "cm")),
                          ...
) {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
        axis.line = axis.line,
        axis.line.x = axis.line.x,
        axis.line.y = axis.line.y,
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.justification = legend.justification,
        legend.title = legend.title,
        legend.box = legend.box,
        panel.background = panel.background,
        panel.border = panel.border,
        panel.grid = panel.grid,
        plot.background = plot.background,
        strip.background = strip.background,
        strip.placement ="outside",
        strip.text = strip.text,
        strip.text.y = strip.text.y,
        ...
    )
}

#' \name{theme_elegant_grey}
#' \alias{theme_elegant_grey}
#' \title{Custom ggplot theme}
#' \description{A customized ggplot theme}
#' \usage{
#' theme_elegant_grey(base_size, base_family, legend.position,
#' legend.title, legend.justification, legend.direction, legend.box, ...)
#' }
#' \arguments{
#' \item{base_size}{The base font size.}
#' \item{base_family}{The base font family.}
#' \item{legend.position}{The position of the legend.}
#' \item{legend.title}{The title of the legend.}
#' \item{legend.justification}{The justification of the legend.}
#' \item{legend.direction}{The direction of the legend.}
#' \item{legend.box}{The box of the legend.}
#' \item{...}{Other arguments passed to `theme`.}
#' }
#' \details{
#' A customized ggplot theme.
#' }
#' \value{A `ggplot` theme.}

theme_elegant_grey <- function(base_size = 11,
                               base_family = "Helvetica",
                               legend.position = "bottom",
                               legend.title = element_text(hjust = 0.5),
                               legend.justification = "center",
                               legend.direction = "horizontal",
                               legend.box = "vertical",
                               ...
) {
  theme(
    text = element_text(family = base_family, size = base_size),
    legend.position = legend.position,
    legend.title = legend.title,
    legend.justification = legend.justification,
    legend.direction = legend.direction,
    legend.box = legend.box,
    ...)
}

# DOCUMENTATION FOR C++ CODE
#-------------------------------------------------------------------------
#' \name{getIndicesCpp}
#' \alias{getIndicesCpp}
#' \title{Retrieve the pair of indices for linear interpolation.}
#' \usage{getIndicesCpp(...)}
#' \description{Retrieve the indices for the interpolation of the
#' LMS parameters
#' }
#' \details{
#' The references for the raw measurements only are given for discrete
#' months values, such that for ages in-between, the L, M and S parameters
#' must be interpolated. For simplicity, the values are linearly interpolated.
#' For this reason, the indices for the lower and upper bound in the proximity
#' of the the values of interest must be found by this function.
#' }
#' \arguments{
#' \item{...}{Arguments.}
#' }
#' \value{
#'  The relevant indices for interpolation.
#' }

getIndicesCpp <- function(...) {
  return(getIndicesCpp(...))
}

#' \name{getProportionCpp}
#' \alias{getProportionCpp}
#' \title{Calculate the proportion factor for precise interpolation.}
#' \usage{getProportionCpp(...)}
#' \description{
#'  Get the relative proportion, which is necessary for interpolation.
#' }
#' \details{
#' The references for the raw measurements only are given for discrete
#' months values, such that for ages in-between, the L, M and S parameters
#' must be interpolated. For simplicity, the values are linearly interpolated.
#' For this reason, the indices for the lower and upper bound
#' in the proximity of the the values of interest must be found by this
#' function.
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \value{
#' The proportion of the difference between the lower / upper
#' reference age bound and the difference of the the age of the child to the
#' upper reference age bound.
#' }
#' \seealso{See also \code{\link{getIndicesCpp}}}

getProportionCpp <- function(...) {
  return(getProportionCpp(...))
}

#' \name{calculateZScoreCpp}
#' \alias{calculateZScoreCpp}
#' \title{z-score transformation}
#' \usage{calculateZScoreCpp(...)}
#' \description{Transform raw biometrical measurements to z-scores}
#' \details{
#' Given the LMS parameters and raw measurements for weight, height or BMI,
#' the given values are transformed to z-values, based on the
#' LMS method of Cole (1991), which is commonly applied for
#' growth curve analysis.
#' }
#' \arguments{
#'   \item{...}{Arguments.}
#' }
#' \value{`vector` of transformed values.}
#' \references{Cole, TJ,
#' "The LMS method for constructing normalized growth standards",
#' European journal of clinical nutrition 44, 1 (1990), pp. 45-60.}

calculateZScoreCpp <- function(...) {
  return(calculateZScoreCpp(...))
}

#' \name{calculateZscoresCpp}
#' \alias{calculateZscoresCpp}
#' \title{z-score transformation}
#' \usage{calculateZscoresCpp(...)}
#' \description{Transform raw biometrical measurements to z-scores}
#' \details{
#' Given the LMS parameters and raw measurements for weight, height or BMI,
#' the given values are transformed to z-values, based on the
#' LMS method of Cole (1991), which is commonly applied for
#' growth curve analysis.
#' }
#' \arguments{
#'   \item{...}{Arguments.}
#' }
#' \value{`vector` of transformed values.}
#' \references{Cole, TJ,
#' "The LMS method for constructing normalized growth standards",
#' European journal of clinical nutrition 44, 1 (1990), pp. 45-60.}
calculateZscoresCpp <- function(...) {
  return(calculateZscoresCpp(...))
}

#' \name{findAgeIndicesCpp}
#' \alias{findAgeIndicesCpp}
#' \title{Find the indices for the interpolation of the LMS parameters.}
#' \usage{findAgeIndicesCpp(...)}
#' \description{Find the indices for the interpolation of the LMS parameters.}
#' \details{
#' The references for the raw measurements only are given for discrete
#' months or years values, such that for ages in-between, the L, M and S parameters
#' must be interpolated. For simplicity, the values are linearly interpolated.
#' For this reason, the indices for the lower and upper bound
#' in the proximity of the the values of interest must be found by this
#' function.
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \value{
#' The relevant indices for interpolation.
#' }

findAgeIndicesCpp <- function(...) {
  return(findAgeIndicesCpp(...))
}

#' \name{calculateCentilesCpp}
#' \alias{calculateCentilesCpp}
#' \title{Calculation of centiles.}
#' \usage{calculateCentilesCpp(...)}
#' \description{Calculates the centiles for a vector of z-values.}
#' \details{
#' This function converts a vector of z-values to centiles based on the LMS method of Cole (1991).
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \value{`vector` of transformed values.}
#' \seealso{See also \code{\link{calculateCentileCpp}}, \code{\link{calculateZScoreCpp}}}
calculateCentilesCpp <- function(...) {
  return(calculateCentilesCpp(...))
}

#' \name{calculateCentileCpp}
#' \alias{calculateCentileCpp}
#' \title{Single centile calculation.}
#' \usage{calculateCentileCpp(...)}
#' \description{Calculate a single centile based on given LMS parameters.}
#' \details{
#' This function deterines the centile based on the given LMS parameters using
#' the LMS method of Cole (1991), which is commonly applied for growth curve analysis
#' based on the formula `centile = median * (1 + L * S * z)^(1 / L)`.
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \value{A single centile value.}
#' \seealso{See also \code{\link{calculateCentilesCpp}}, \code{\link{calculateZScoreCpp}}}
calculateCentileCpp <- function(...) {
  return(calculateCentileCpp(...))
}

#' \name{assignMedianBmiCpp}
#' \alias{assignMedianBmiCpp}
#' \title{Assign the median BMI to the children.}
#' \usage{assignMedianBmiCpp(...)}
#' \description{Assign the median BMI to the children.}
#' \details{
#' Assign the median BMI to the children.
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \value{`vector` of transformed values.}
#' \seealso{See also \code{\link{findAgeIndicesCpp}}}

assignMedianBmiCpp <- function(...) {
  return(assignMedianBmiCpp(...))
}

#' \name{interpolateCpp}
#' \alias{interpolateCpp}
#' \title{Linear interpolation between an upper and lower value.}
#' \usage{interpolateCpp(...)}
#' \description{Perform linear interpolation based on the lower / upper
#' reference age bound.}
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \details{
#' The basic functionality to perform linear interpolation using the proportion
#' of the upper and lower age bound,
#' relevant for the respective age of a child.
#' }
#' \value{A 'numeric' value, representing the linear interpolation.}
#' \seealso{See also \code{\link{getIndicesCpp}}, \code{\link{getProportionCpp}},
#' \code{\link{interpolateToZscoreCpp}}}

interpolateCpp <- function(...) {
  return(interpolateCpp(...))
}

#' \name{interpolateToZscoreCpp}
#' \alias{interpolateToZscoreCpp}
#' \title{Linear interpolation of LMS parameters.}
#' \usage{interpolateToZscoreCpp(...)}
#' \description{
#'  Interpolation of a vector of raw values of either weight, height or BMI to z-scores.
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \details{
#' For getting z-scores for values for ages, not covered by references tables,
#' linear transformation for ages in close proximity are used.
#' }
#' \value{
#' A 'vector' of z-score transforms.
#' }
#' \seealso{See also \code{\link{getIndicesCpp}}, \code{\link{getProportionCpp}}, \code{\link{interpolateCpp}}.}

interpolateToZscoreCpp <- function(...) {
  return(interpolateToZscoreCpp(...))
}

#' \name{getAgeCohortsCpp}
#' \alias{getAgeCohortsCpp}
#' \title{Assign age bounds to observations.}
#' \usage{getAgeCohortsCpp(...)}
#' \description{
#' Get the age bounds for each observation based on the age of the subject.
#' }
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \details{
#' For each observation the age is compared to the given age bounds, such that the given age range
#' is assigned to the respective observation, if the age at the time of measurement lies in the
#' given range. Note, that if the age cannot be assigned to any range, the assigned value is zero.
#' }
#' \value{
#' A list of with upper and lower age bounds.
#' }
#' \seealso{See also \code{\link{getIndicesCpp}}, \code{\link{getProportionCpp}}, \code{\link{interpolateCpp}}.}

getAgeCohortsCpp <- function(...) {
  return(getAgeCohortsCpp(...))
}

#' \name{getUniqueIdsCpp}
#' \alias{getUniqueIdsCpp}
#' \title{Retrieve unique IDs.}
#' \usage{getUniqueIdsCpp(...)}
#' \description{Return the unique IDs for a all observations.}
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \details{
#' All unique subject IDs are extracted.
#' }
#' \value{
#' A `vector` of unique subject names.
#' }
#' \seealso{See also \code{\link{getIndicesCpp}}, \code{\link{getProportionCpp}}, \code{\link{interpolateCpp}}.}

getUniqueIdsCpp <- function(...) {
  return(getUniqueIdsCpp(...))
}

# Find the age values for the reference table
.mapAgeValuesCpp <- function(...) {
  return(mapAgeValuesCpp(...))
}

#' \name{findIndicesCpp}
#' \alias{findIndicesCpp}
#' \title{Find the indices for the age values in the reference table.}
#' \usage{findIndicesCpp(...)}
#' \description{Find the indices for the age values in the reference table.}
#' \arguments{
#'  \item{...}{Arguments.}
#' }
#' \details{
#' Find the indices for the age values in the reference table.
#' }
#' \value{
#' A `vector` of indices.
#' }
#' \seealso{See also \code{\link{getIndicesCpp}}.}

# Find the indices for the age values in the reference table
findIndicesCpp <- function(...) {
  return(findIndicesCpp(...))
}

# Check the reference name
.check_ref_name <- function(ref_name) {
  if (!any(tolower(ref_name) %in% tolower(get_reference_names()))) {
    stop("No reference with name '", ref_name, "' found.")
  }
}

# Check the measure
.check_measure <- function(measure) {
  if (!all(tolower(measure) %in% tolower(c("Weight", "Height", "BMI")))) {
    stop ("measure must be one of 'Weight', 'Height' or 'BMI'")
  }
}

# Helper function to allocate variables used in a data.table
# Takes a 'character' or 'lang' vector of variable names and assigns NULL
.allocate <- function(x) {

  # This also creates the variables based on their
  # names into the scope of the parent, i.e. usually the function,
  # .allocate was called from
  for (v in x) {
    assign(v, NULL, envir = parent.frame())
  }

  return(x)
}

# R packages should contain .rda / .RData files
# but these are inflexible, i.e. you have to use
# the exact name of the file. The function below
# makes sure that the right file is loaded
.load_rda <- function(rda_file) {
  rda_name <- str2lang(load(rda_file))
  return(eval(rda_name))
}

# Function to map the reference tables to the correct measures
.map_ref <- function(ref_name, measure) {
  # This retrieves the absolute path of the .rda files, which are located in inst/extdata
  extdata <- system.file("extdata", package = getPackageName())

  filenms <- if (measure != "all") {
    paste0(extdata, "/", tolower(ref_name), "_", tolower(measure), ".rda")
  } else {
    paste0(extdata, "/", tolower(ref_name), ".rda")
  }
  filenms <- filenms[file.exists(filenms)]

  if (length(filenms) == 0) {
    stop("Reference file not found")
  } else if (length(filenms) == 1) {
    return(.load_rda(filenms))
  } else {
    ref <- .load_rda(filenms[1])
    for (fn in filenms[-1]) ref <- merge(ref, .load_rda(fn), all = TRUE)
  }

  return(ref)
}

# TODO: This should be checked
# Primitive interpolation function
.interpl <- function(x, y = NULL, xout, interpol_method = "linear", spline_method = "natural", ...) {
  switch(interpol_method,
         "linear" = approx(x = x, y = y, xout = xout, ...)$y,
         "spline" = spline(x = x, y = y, xout = xout, method = spline_method, ...)$y
  )
}

###  place startup loads here
.onLoad <- function(libname, pkgname) {
  # to show a startup message
  # example on how to initialize_data a Tcl package
  # tcltk::.Tcl(paste("lappend auto_path",file.path(system.file(package="Rpkg"),"pantcl", "lib")))
  # tcltk::.Tcl("package require tclfilters")
  # tools::vignetteEngine("pantcl",package=pkgname,weave=pantcl,tangle=ptangle)
  # Expose the interpolation module to R
  Rcpp::loadModule("interpolation_module", TRUE)
}
