---
title: "evaluating_gpbart"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{evaluating_gpbart}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




### Setting a example

This a vignette to explain how to run a simple example of the model, setting its own prior and its hyperparameters. To start we going to use the `friedman` example as the dataset to be used.


```r
library(testgpbart)

# Setting simualation parameters
n <- 100
seed <- 42
sd <- 0.1

# Loading the data
fried_data <- sim_friedman(n = n,seed = seed,sd = sd)

# Setting a cross-validation object
cv_obj <- k_fold(data = fried_data,dependent_variable = "y",
                 k_partitions = 5,seed = seed)

# Selecting one fold
fold <- 1
x_train <- cv_obj[[fold]]$x_train
y_train <- cv_obj[[fold]]$y_train
x_test <- cv_obj[[fold]]$x_test
y_test <- cv_obj[[fold]]$y_test

```

Once specified all the data setting the next necessary specification is with respect the model settings. Here the most important thing that we want to evaluate is the 
prior for $\phi_{tp}$. To specify that there are two main arguments in the function: `proposal_phi_` and `prior_phi_` that are going to specify the proposal and the prior for a MH sampling in that algorithm. For the `proposal_phi_` we need to enter a list with the `proposal_mode` and a possible `grid` to be used depending on the defined proposal. Examples for each case are explained below

- **Sliding-window proposal**: for this proposal the proposal for a new $\phi^{*}_{tp}$ came from a sample that follows $$\phi^{*}_{tp} \sim U\left(3/4\phi_{tp},4/3\phi_{tp}\right)$$. Therefore there is no grid to be defined and the function argument should be given by 

```r
proposal_phi_ = list(proposal_mode = "sliding_window", grid = NULL)
```

- **Discrete grid**: for this proposal be given by a discrete uniform distribution with a range of $\phi^{*}_{tp}$ values. Therefore, is necessary to enter a vector with all values that $\phi_{tp}$ in that grid. The example below demonstrate one option of possible grid to be used


```r
proposal_phi_ = list(proposal_mode = "discrete_grid",
                     grid = c(seq(0.1,20,length.out = 20),100))
```

Lastly, the default argument of the function takes uses


```r
proposal_phi_ = list(proposal_mode = "discrete_grid",
                     grid = NULL)
```

where a default grid (showed below) is used as the standard one.


```r
phi_proposal_ <- sample(c(0.1,seq(0,10,by=0.5)[-1],seq(10,20,by = 1),
                         25,30,50,75,100,125),
                       size = 1)
```


For the prior definition we have that it will gonna follow a mixture of gammas $$ \pi_{1}G(\alpha = \alpha_{1}, \beta = \beta_{1})+\pi_{2}\left(\alpha = \alpha_{1}, \beta = \beta_{1} \right)$$, where $\pi_{i}$ corresponds to the probability of each gamma, and $\alpha_{i}$ and $\beta_{i}$ are the shape and the rate parameters of each one respectively. Therefore we need to specify all parameters of it through the list


```r
phi_prior_ <- list(type = "gamma_mixture",
                   prob_1 = 0.5,shape_1 = 1,rate_1 = 1,
                   prob_2 = 0.5,shape_2 = 1000,rate_2 = 10)
```

The list above would refer to a prior given by $$ 0.5G(\alpha = 1, \beta = 1)+\pi_{2}\left(\alpha = 1000, \beta = 10 \right)$$


If null arguments are chosen for the `phi_prior_` the default prior is used follows
$$ 0.3 \times G(\alpha = 1.5, \beta = 0.2)+0.7 \times G\left(\alpha = 10000, \beta = 100 \right)$$
The code for it would be


```r
proposal_phi_ = list(type = NULL,
                     prob_1 = NULL, shape_1 = NULL,rate_1 = NULL,
                     prob_2 = NULL,  shape_2 = NULL, rate_2 = NULL)
```

### Running the model

Finally to run the model we would have:


```r
gp_bart_mod <- gp_bart(x_train = x_train,
                       y_train = c(y_train),
                       x_test = x_test,
                       n_tree = 10,bart_boolean = TRUE,
                       gp_variables_ = colnames(x_train), # Selecting all var.
                       rotation_variables_ = colnames(x_train)) 
```

### Real data benchmarking

To verify over other real data benchmarking use the other data from the package


```r
# All other datasets
auckland
baltimore
boston
columbus
ny_PCTOWNHOME
ny_PROPCAS
sponge
swmud
petrel
precip
```