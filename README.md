# SSPTools#

SSPTools is a MATLAB toolkit for investigating the Strong Stability Preserving (SSP) properties of numerical methods that
approximate solutions to conservative hyperbolic partial differential equations. 

More information about SSP methods and current researh projects can be found at the http://www.sspsite.org

Please check out [SSP_Tools](https://github.com/DanielHiggs/SSP_Tools),the first version written by Dan Higgs.

## Integrator
* Runge-Kutta Method (using the Butcher Tableau)
	- ERK (Explicit Runge-Kutta)
	- DIRK (Diagonally Implicit Runge-Kutta)

## Discretizers
* Finite Difference Method
	- Specify stencil direction 
		1. CD - Central Difference
		2. FD - Forward Difference
		3. BD - Backward Difference
	- [Now allows for different order of accuracy, for either directions](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
* NDG (Nodal Discontinuous Galerkin)
	- Only tested for:
		1. Linear Advection
		2. Burgers
* WENO ( Weighted Essentially
Non-Oscillatory Schemes)
	- Only WENO5 is working for:
		1. Advection
		2. Burgers

## Tests
* Convergence
* SSP

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
	- please write unit tests and make sure that all previously written unit tests are running and `passing!!`
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact

### Authors ###
The following people have contributed code to SSPTools 

 - [Sidafa Conde](http://hilbert.math.umassd.edu/~sconde/): Principal author and maintainer
