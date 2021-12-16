# CVXSynth

Implementations of the original synthetic control (Abadie, Diamond, Hainmueller 2010, 2015) and elastic net synthetic control (Doudchenko and Imbens 2016).

+ [vignette](https://apoorvalal.github.io/posts/09122021_ElasticNetSyntheticControl.html)
+ [notes](https://apoorvalal.github.io/presentations/pdf/ImbensDoudchenko.pdf)

## Installation

```{r}
library(remotes)
install_github("apoorvalal/CVXSynth")
```

## Optimiser

Uses the domain-specific language [CVXR](https://cvxr.rbind.io/) with the
commercial solver `MOSEK` to solve for weights. `MOSEK` is free for academics
and can be installed using installation instructions on their website. Else, an
[alternate
solver](https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/) can be
passed as the `solv` argument for both `sc_solve` and `en_sc_solve` functions.
