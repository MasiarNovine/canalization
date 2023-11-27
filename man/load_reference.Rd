\name{load_reference} 
\alias{load_reference} 
\title{Load reference tables.} 
\description{Using reference tables} 
\arguments{ 
 \item{ref_name}{Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".} 
 \item{measure}{A vector of anthrometric measure to be used. Default: c("Weight", "Height", "BMI")} 
 \item{na_omit}{Should NAs be omitted? Logical. Default: FALSE.} 
} 
\details{ 
The World Health Organization (WHO) provides an international composite reference based on 
several national references from different countries and continents. 
`KiGGS` is a reference based on the German nation-wide KiGGS published by the Robert-Koch Institute in Germany. 
`Kromeyer-Hauschild` is a another German national reference. Note, that currently, only WHO and KiGGS provide 
LMS parameters for all measures, while `IOTF` denotes only LMS parameters for BMI and `Kromeyer-Hauschild` does 
not include L and S parameters for weight. 
} 
\value{A reference as a `data.table` object.} 
\references{ 
TJ Cole, "The LMS method for constructing normalized growth standards", European journal of clinical nutrition 44, 1 (1990), pp. 45-60, 
Stolzenberg, H., H. Kahl, and K. E. Bergmann (May 2007). “Körpermaße bei Kindern und Jugendlichen in Deutschland”. In: Bundesgesundheitsblatt - Gesundheitsforschung - Gesundheitsschutz 50.5, pp. 659–669. issn: 1437-1588. url: https://doi.org/10.1007/s00103-007-0227-5, 
Kromeyer-Hauschild, K. et al. (Aug. 2001). “Perzentile für den Body-mass-Index für das Kindes- und Jugendalter unter Heranziehung verschiedener deutscher Stichproben”. In: Monatsschrift Kinderheilkunde 149.8, pp. 807–818. issn: 1433-0474. url: https://doi.org/10.1007/s001120170107. 
} 
