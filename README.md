# spritzer
Search for periodic signals in noisy time series using Fourier methods

spritzer is a set of R functions developed for detecting strictly or nearly periodic signals in time series dominated by red noise. It is based on code developed and used in

**[A Bayesian test for periodic signals in red noise](http://adsabs.harvard.edu/abs/2009arXiv0910.2706V)** 
Vaughan S., 2010, MNRAS, 402, 307

The name SPRITZER is a blend of the words Strictly Periodic Tester.

## Installation

spritzer is an R package, but is still in development. To set up from GitHub first install Hadley Wickham's devtools.
```
install.packages("devtools")
```
Now you can install straight from GitHub:
```
devtools::install_github("svdataman/spritzer")
```
It requires the mnormt package. If you don't already have this, you will also need to install it locally:
```
install.packages("mnormt")
require(spritzer)
```
and you're good to go.

## Basic usage

Given a data file as input, containing the time series data in columns (e.g. time, value, error [optional]) of a plain text file

```
result <- spritz("data/mrk766.txt")
```

## Assumptions

Spritzer works best if the input time series is regularly sampled, with no gaps, and is a realisation of a red noise process (with a smooth, steep power spectrum), which might contain an additional strictly periodic component. See my 2010 paper for more details.

## To do

Lots to do before this is finished. 

* replace the MCMC engine with one from tonic (chain generation and diagnostic plots)
* strip out the file loading, assume user has data in memory as array (or data.frame)
* better format of output list
* strip out the interactive element (user response to questions) and replace with function arguments on inputs.
* replace simulation in Fourier space with full simulations in time. 

## License

MIT
