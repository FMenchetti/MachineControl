
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MachineControl

<!-- badges: start -->
<!-- badges: end -->

## What is MachineControl?

MachineControl is an R package that enables the estimation of both the average and heterogeneous causal effects of a simultaneous policy or shock (also referred to as ‘treatment’ or ‘intervention’) that directly or indirectly impacts all the units in a panel dataset. This presents a challenging causal inference scenario, as there are no unaffected units available for constructing a valid comparison group. In such situations, traditional methods like Difference-in-Differences (DiD) or the Synthetic Control Method (SCM) cannot be applied, as they rely on control units to estimate the counterfactual outcome.

## How does the package work?

MachineControl operates by harnessing Machine Learning (ML) algorithms to understand the temporal dynamics and associations between the units and covariates in the pre-intervention panel dataset. The process involves a horse-race competition among various ML methods: all the ML
methods selected by the user are tuned via a panel cross-validation
procedure in order to track as best as possible the pre-intervention trend of the outcome variable. Then, the procedure selects the best-performing algorithm
based on a given performance metric, such as the Root Mean Squared Error
(RMSE).

Once the best-performing model is selected, it is used to forecast the counterfactual outcome for each unit—essentially, what would have occurred in the post-intervention period had the treatment been absent. The unit-level difference between the observed outcome and the counterfactual outcome represents the estimated individual treatment effect, which can then be aggregated to get the Average Treatment Effect (ATE).
MachineControl can also be employed to learn about treatment effect heterogeneity by estimating data-driven Conditional Average Treatment Effects (CATEs). 

Inference on the estimated effects (ATE and CATE) is conducted using a block-bootstrapping approach.

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

Cerqua A., & Letta M., & Menchetti F. 2023. The Machine
Learning Control Method for Counterfactual Forecasting. Available as SRRN Working Paper at: http://dx.doi.org/10.2139/ssrn.4315389
