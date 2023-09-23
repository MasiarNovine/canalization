# canalization

This R package contains the source code used in the context of the Master's thesis at
the University of Potsdam.

## Building the package

The package can be build by using the `make` utility:

```
make build
```

This runs a small R function (provided by Detlef Groth), which extracts the documentations from the
source file and creates `Rd` files for each documented function. After that, it builds and checks the tar-ball.

If building has worked, the package can be installed by the common `R CMD INSTALL <tar-ball-name>`.
