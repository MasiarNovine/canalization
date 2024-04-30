# canalization

This repository contains the source code of the R package used for the analysis part of my Master's thesis at
the University of Potsdam. 

The analysis was concerned with the investigation of the relationship of BMI volatility and the risk for obesity based on longitudinal data provided by CrescNet (University of Leipzig). The data consisted of randomly stratified sampled German children and adolescents from the years 2000 to 2022.     

## Building the package

The package can be build by using the `make` utility:

```
make build
```

This runs a small R function (provided by Detlef Groth), which extracts the documentations from the
source file and creates `Rd` files for each documented function. After that, it builds and checks the tar-ball.

If building has worked, the package can be installed by the common `R CMD INSTALL <tar-ball-name>`.
