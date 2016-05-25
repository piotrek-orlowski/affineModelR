# affineModelR

R package for working with multifactor stochastic volatility models, as in Duffie, Pan and Singleton (2000) 'Transform Analysis and Asset Pricing for Affine Jump-Diffusions'

# Installation

The package has been checked to compile from source on Linux and Windows. The `-fpermissive` compilation flag is necessary for using a compiled function in the ODE solver (`deSolve` package).

To install the package, run:
```
library(devtools)
install_github(repo = "piotrek-orlowski/affineModelR")
```

# Use cases

This package allows for calculating the values of the characteristic function in Affine Jump Diffusion models with an arbitrary number of factors and flexible jump specifications. The CF can be evaluated both under the statistical probability measure *P*, and the risk-neutral (pricing) measure *Q*. Knowledge of the characteristic function is equivalent to the knowledge of the distribution of a random variable.

The *P* CF can be used for calculating conditional moments of stock returns and volatility factors, for example for purposes of parameter estimation, filtering. The *Q* CF can be used for calculating prices of derivatives with various transform methods (the package https://github.com/piotrek-orlowski/transformOptionPricer includes a basic functionality for vanilla options).

The package also contains a fast simulating engine: it is easy to draw large samples from AJDMs at arbitrary frequencies. This feature can be useful in advanced courses in Derivatives, or for developing tools in High-Frequency Econometrics.

# Examples

The file https://github.com/piotrek-orlowski/affineModelR/blob/master/inst/tests/examples-CF.R contains examples of calculations that the package allows for.

# External headers

The package makes some of its functions available via headers in `inst/include`. These functions include (1) c++ evaluators of CFs given coefficients returned from `affineCF` or `affineCFderivs` in file `affineCF.h` and (2) typedefs for providing users jump transforms to the affine ODE solver in `jumpTransforms.h`. Users who would like to use these functionalities, have to build their packages against `affineModelR` by providing appropriate flags in the `Makevars` files. The https://github.com/piotrek-orlowski/divergenceModelR package installation process explains this in more detail.

# Your own jump specification

The package now supports two types of co-jump specifications. First, it is assumed that only the first stochastic volatility factor jumps alongisde the stock price. In both cases the jumps in volatility are exponentially distributed. The location parameter of the jump in the asset price depends on the realisation of the volatility jump. Conditional on the latter, the asset price jump follows either a Gaussian, or a Laplace distribution.

If you want to add your own jump specification, you have to take care of *jump simulation* and *moment calculation* separately. In order to add a jump generator or a jump transform, `#include` the header `affineModelR.h` in your code. In order to add a jump specification to CF calculation, you have two choices. If you only want to calculate the CF values (`affineCF`), you have to write a function that accepts `arma::cx_colvec` (argument) and `Rcpp::List` (parameters) and returns `std::complex<double>`, then pass it as an `Xptr` to `affineCF`. If you would like to use `affineCFderivs`, you have to provide similar functions which provide derivatives (gradient, hessian, gradient-of-hessian, respectively) of your jump transform. Consult `jumpTransform.h` and `jumpTransform.cpp` for details. 

Your jump transform has to accept an `arma::vec` of parameters and return an `arma::vec` of values, one draw per call. Then pass the function as an `Xptr` to `affineSimulate`. Simulation also requires providing the base jump transform function pointer.

This admittedly complex mechanism allows for flexibility in re-specifying the models while retaining the speed of the `C++`-based ODE solver.

# Authors

The package has been developed by Andras Sali (https://github.com/andrewsali) and Piotr Orłowski.
