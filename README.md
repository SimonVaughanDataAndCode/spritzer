# spritzer
Search for periodic signals in noisy time series using Fourier methods

spritzer is a set of R functions developed for detecting strictly or nearly periodic signals in time series dominated by red noise. It is based on code developed and used in

**A Bayesian test for periodic signals in red noise** 
Vaughan S., 2010, MNRAS, 402, 307
[(http://adsabs.harvard.edu/abs/2009arXiv0910.2706V)]

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

## Usage

Given a data file as input, containing the time series data in columns (e.g. time, value) of a plain text file

```
result <- bayes("data/mrk766.txt")
```

