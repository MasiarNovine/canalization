\name{load_ref} 
\alias{load_ref} 
\title{Load references.} 
\description{Loads LMS references.} 
\usage{load_ref(ref_name = "Kromeyer-Hauschild", measure = c("Weight", "Height", "BMI", "all"))} 
\arguments{ 
  \item{ref_name}{Name of the reference table.} 
  \item{measure}{Vector of measures.} 
} 
\details{ 
Loads reference tables for measures `Weight`, `Height` and `BMI`. 
The reference values used are based on the GAMLSS fits provided by the 
\code{childsds} package (author: M.Vogel), such that the parameter names 
are nu (L), mu (M) and sigma (S). 
} 
\value{`data.table` with nu, mu and sigma parameters.} 
\references{ 
TJ Cole, "The LMS method for constructing normalized growth standards", European journal of clinical nutrition 44, 1 (1990), pp. 45-60, 
Stolzenberg, H., H. Kahl, and K. E. Bergmann (May 2007). “Körpermaße bei Kindern und Jugendlichen in Deutschland”. In: Bundesgesundheitsblatt - Gesundheitsforschung - Gesundheitsschutz 50.5, pp. 659–669. issn: 1437-1588. url: https://doi.org/10.1007/s00103-007-0227-5, 
Kromeyer-Hauschild, K. et al. (Aug. 2001). “Perzentile für den Body-mass-Index für das Kindes- und Jugendalter unter Heranziehung verschiedener deutscher Stichproben”. In: Monatsschrift Kinderheilkunde 149.8, pp. 807–818. issn: 1433-0474. url: https://doi.org/10.1007/s001120170107. 
Rigby, R.A. and D.M. Stasinopoulos (2005). “Generalized additive models for location, scale and shape.” In: Applied Statistics 54, pp. 507–554. 
} 
