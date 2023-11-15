
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MachineControl

<!-- badges: start -->
<!-- badges: end -->

## What is MachineControl?

MachineControl is an R package that can be used to estimate both the
average and the heterogeneous causal effect of a simultaneous policy or shock (also called
‘treatment’ or ‘intervention’) that affects, directly or indirectly, all the units
of a panel dataset. This is a challenging causal inference setting as
there are no unaffected units that can be employed to build a valid comparison group. In such settings, traditional methods such as
Difference-in-Differences (DiD) or the Synthetic Control Method (SCM) cannot be applied, as they
require control units to estimate the counterfactual outcome in the
absence of the treatment.

## How does the package work?

MachineControl works by using a pool of Machine Learning (ML) algorithms to
learn the temporal dynamics and the associations between the units and
the covariates in the pre-intervention panel dataset. It involves a
horse-race between different ML methods: all the ML
methods selected by the user are tuned via a panel cross-validation
procedure in order to give the best predictions of the outcome before
the intervention. Then, the procedure selects the best-performing method
based on a given performance metric, such as the Mean Squared Error
(MSE).

Once the best-performing model is selected, it is used to forecast the counterfactual outcome for each unit, i.e., what
would have happened in the post-intervention period absent the treatment. The difference between the observed
outcome and the counterfactual outcome is the estimated Average
Treatment Effect (ATE) of the policy. MachineControl can also be used to
estimate the heterogeneous effect of the treatment. For example, while
estimating the impact of COVID-19 on mortality, regions having a lower
number of hospital beds per 1000 people or with a higher proportion of
elderly population could be more affected than others. This is usually
referred as Conditional Average Treatment Effect (CATE), as we are
interested in averaging the effects on subsets of units sharing similar
characteristics.

Inference on the estimated effects (ATE and CATE) is performed by
block-bootstrapping approach.

## Installation

You can install the development version of MachineControl from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("FMenchetti/MachineControl")
```

<!--
## Example
&#10;This is a basic example which shows you how to solve a common problem:
&#10;
```r
library(MachineControl)
#> Loading required package: caret
#> Loading required package: ggplot2
#> Loading required package: lattice
## basic example code
```
-->

## Further readings

Cerqua A., & Letta M., & Menchetti F. The Machine
Learning Control Method for counterfactual forecasting. Available as SRRN Working Paper at: http://dx.doi.org/10.2139/ssrn.4315389
