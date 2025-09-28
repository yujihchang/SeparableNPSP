**Separable Effects for (Semi)Competing Risks**

This R package implements two approaches for estimating separable direct and indirect effects in mediation analysis under competing risks and semicompeting risks:

- **Nonparametric method (NP):**  
  A model-free estimator based on multistate hazard decomposition and Nelson–Aalen cumulative hazard increments.  
  *Advantage:* robust without model assumptions.  
  *Limitation:* less efficient when adjusting for confounders.

- **Semiparametric method (SP):**  
  A Cox proportional hazards–based estimator that flexibly adjusts for baseline confounders by including them as covariates.  
  *Advantage:* efficient with covariate adjustment.  
  *Limitation:* relies on proportional hazards assumption.

Both methods are motivated by the paper:  

> Yu, J.-C. and Huang, Y.-T. (2025).  
> *Separable Effects of Semicompeting Risks: The Effects of Hepatitis B on Liver Cancer via Liver Cirrhosis.*  
> *Statistics in Medicine.* doi: [10.1002/sim.70178](https://doi.org/10.1002/sim.70178)

---

## Installation

```r
# From GitHub (example)
devtools::install_github("your-org/Separable_NPSP")
