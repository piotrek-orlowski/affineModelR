# affineModelR

R package for working with multifactor stochastic volatility models, as in Duffie, Pan and Singleton (2000) 'Transform Analysis and Asset Pricing for Affine Jump-Diffusions'

# Installation

The package has been checked to compile from source on Linux and Windows. The `-fpermissive` compilation flag is necessary for using a compiled function in the ODE solver (`deSolve` package).

# Use cases

This package allows for calculating the values of the characteristic function in Affine Jump Diffusion models with an arbitrary number of factors and flexible jump specifications. The CF can be evaluated both under the statistical probability measure *P*, and the risk-neutral (pricing) measure *Q*. Knowledge of the characteristic function is equivalent to the knowledge of the distribution of a random variable.

The *P* CF can be used for calculating moments of stock prices and volatility factors, for example for purposes of parameter estimation. The *Q* CF can be used for calculating prices of derivatives with various transform methods (the package https://github.com/piotrek-orlowski/transformOptionPricer includes a basic functionality for vanilla options).

The package also contains a fast simulating engine: it is easy to draw large samples from AJDMs at arbitrary frequencies. This feature can be useful in advanced courses in Derivatives, or for developing tools in High-Frequency Econometrics.

# External headers

The package makes some of its functions available via headers in `inst/include`. The end user has to indicate to the compiler how to link against the package.

# Authors

The package has been developed by Andras Sali (https://github.com/andrewsali) and Piotr Or{\l}owski.
