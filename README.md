# canalization

This R package contains the source code used for the analysis part of my Master's thesis at
the University of Potsdam. The analysis was concerned with the investigation of the association of BMI volatility and the risk for obesity, using a dataset of German children and adolescents from the years 2000 to 2022 provided by CrescNet, a research network associated with the University of Leipzig.    

## Building the package

The package can be build by using the `make` utility:

```
make build
```

This runs a small R function (provided by Detlef Groth), which extracts the documentations from the
source file and creates `Rd` files for each documented function. After that, it builds and checks the tar-ball.

If building has worked, the package can be installed by the common `R CMD INSTALL <tar-ball-name>`.
