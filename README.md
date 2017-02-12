# spritzer
Search for periodic signals in noisy time series using Fourier methods

spritzer is a set of R functions developed for detecting strictly or nearly periodic signals in time series dominated by red noise. It is based on code developed and used in

**[A Bayesian test for periodic signals in red noise](http://adsabs.harvard.edu/abs/2009arXiv0910.2706V)** 
Vaughan S., 2010, MNRAS, 402, 307

The name SPRITZER is a blend of the words Strictly Periodic Tester.

## Installation

spritzer is not (yet) an R package, and is still in development. To set up the R functions source the .R files:
```
source("bayes.R")
```

It requires the mnormt package. If you don't already have this, you will also need to install it locally:
```
install.packages("mnormt")
```

## Basic usage

Given a data file as input, containing the time series data in columns (e.g. time, value, error [optional]) of a plain text file

```
result <- bayes("data/mrk766.txt")
```

## Assumptions

Spritzer works best if the input time series is regularly sampled, with no gaps, and is a realisation of a red noise process (with a smooth, steep power spectrum), which might contain an additional strictly periodic component. See my 2010 paper for more details.
