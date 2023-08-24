# Sex conflict evolution
Project to analyze the evolution of sex conflict in new and old genes of Drosophila melanogaster. We will be modeling offspring distributions in non-essential gene knock-down lines. Project consists of egg count data and fitted Poisson or Negative binomial regression using RStan. Posterior estimates of gene knockdown effects in somatic and germline cells, male and female parentals, or cross combinations will be compared to assess intrahost conflict or sex conflict, respectively, and correlated with gene age and function.

## Model goals
We want to identify genes non-essential for development that have apparent conflicting effects on fertility. This conflict would arise as differing fertility effects between somatic and germline expression or between male and female expression of a candidate gene. We also want to identify if adaptive and non-adaptive genes (for males, females, or both) correlate with gene age, chromosome location, and expression pattern. 

## Analysis steps
To address these question, we need to estimate the effects of RNAi knock down experiments on fertility. Using UAS RNAi, we can knock down gene expression in the germline using nos GAL4, knock down gene expression in the soma using Act05 GAL4, and knock down expression for both using Tubp GAL4. 

We want to categorize the mean and shape of the offspring distribution for these lines compared to driver controls. To do so, we'll parameterize statistical models with offspring count data. The formulation of these candidate models will differ to account for various biological assumptions, including the dispersion of the offspring distribution, the source of zero-count data, and sex-specific interactions with those parameters. The best candidate model will be identified using formal model selection to see which formulation offers the most explanatory power for our data. 

From the best model, we will take parameters estimating the effects of our reverse genetics experiments on offspring distribution. These parameters will be used in our analyses relating gene age, chromosome location, and expression pattern to fertility effects as we will claim that they represent the ``true" values of fertility impact by each gene.

# Model structure
## Base one-gene model

The simplest realistic model we can generate for a single gene knock down experiment assumes that the offspring count $M$ follows a Poisson distribution with mean $\mu$ offspring. It is realistic to assume that at least the genetic background, driver, and RNAi knock down affect this mean. We can treat these as coefficients that are multiplied by the base reproduction mean $\beta_0$ to yield the total reproduction mean $\mu$. Since we fit models on a log scale, this looks like

```math
\begin{align}
    \log(\mu) &= \beta_0 + & \text{base reproduction} \\
    & \beta_L + & \text{genetic background line effect} \\
    & \beta_D  + & \text{driver effect} \\
    & \beta_{K;d} & \text{knock down effect for driver $d$}\\
    M&\sim\text{Pois}(\mu)
\end{align}
```

where each $\beta$ term affects the overall mean offspring count $\mu$. The sum of the $\beta$ terms is exponentiated to yield the mean expected number of offspring. 

## Realistic adjustments for simplest model
The base model is reasonable as a single gene model, but the model is so simple compared to our data and later elaborations that it's disingenuous to even attempt fitting it. When confronted with data and compared to other models, we know that it's going to lose. Therefore, the truly simplest model that we'll move forward with is one where the  the knock down effects are not just of course gene-specific and driver-specific but also sex-specific. We will also allow the effect of each driver to be sex-specific and genetic background-specific. Our simplest model then is

```math
\begin{align}
\log(\mu) &= \beta_0 + & \text{base reproduction} \\
    & \beta_L + & \text{genetic background line effect} \\
    & \beta_{D;sl}  + & \text{driver effect for sex $s$ and line $l$} \\
    & \beta_{K;gdsl} & \text{knock down effect for gene $g$, driver $d$, sex $s$, line $l$}\\
    M&\sim\text{Pois}(\mu)
\end{align}
```

and any models that we compare to this simplest model will be more elaborate. 

## Over-dispersion
We can allow for over-dispersion, where an excess of high offspring counts occur, by changing the distribution type and adding an extra parameter. Negative-binomial distributions are well-suited for this kind of data and introduce the dispersion parameter $r$. This looks like the simplest model but offspring count $M$ is distributed differently:

```math
\begin{align}
    \log(\mu)&=\beta_0 +\beta_L + \beta_{D;sl}  + \beta_{K;gdsl} \\
    M&\sim\text{NegBinom}(\mu, r)
\end{align}
```

We can also say that overdispersion $r$ is sex-specific, which would add two parameters, $r_\text{male}$ and $r_\text{female}$, instead of just one. 

## Excess zeroes
There are higher than expected zero counts in many of the crosses, so we need to account for these in some ways. Two common approaches we can compare are zero inflation and a zero hurdle 

### Zero-inflation
Zero-inflation is a mixture distribution that depends on a probability $\theta$. If we say that the offspring distribution model is some function $f(m)$,

```math
\begin{equation}
    P(M=m|\theta,\vec\beta,r)=
    \begin{cases}
        \theta + (1-\theta)\times f_\text{model}(0) & m = 0\\
         (1-\theta)\times f_\text{model}(m) & m > 0\\
    \end{cases}
\end{equation}
```

The probability of observing a zero comes from the model offspring distribution with probability $1-\theta$ and some factor external to the model distribution with probability $\theta$


### Zero hurdle
Zero-hurdle is similar to zero-inflation but the zeroes all come from a factor external to the distribution and we say that the offspring count distribution is truncated to exclude zero. We adjust for this by dividing the offspring distribution by the probability of observing a zero value.

```math
\begin{equation}
    P(M=m|\theta,\vec\beta,r)=
    \begin{cases}
        \theta & m = 0\\
         (1-\theta)\times\frac{f_\text{model}(m)}{f_\text{model}(0)} & m > 0\\
    \end{cases}
\end{equation}
```

In both of these models, we can make multiple values of $\theta_{gsd}$ that depend on the gene of interest, sex, and driver. These specific $\theta$ terms can be included with or without a base level $\theta_0$.

## Experiment covariates
The only observed condition in the vial considered to possibly affect the offspring distribution is early mortality of parentals, which we can quickly add by adding a factor $\beta_P$ to the log mean expectation $log(\mu)$. 

There were observed effects of time of year start time for the experiment. We can simply include various coefficients for either a particular start day where each start day has some value $\beta_t$ or each experiment started at the same month has some value $\beta_T$. These are likely correlated with each other and we can complicate things with time correlation, but these further complications are likely to be difficult to fit.

# Model fitting (too briefly)
We'll parameterize the models above using Markov-chain Monte Carlo (MCMC) in the RStan package. The resulting ``chains" will each be samples from the posterior distributions of our biological parameters based on our priors and likelihood functions. 

Model selection helps to determine the quality of fits made using Bayesian inference, and rank models based on how accurate the predictions they make are. The best models will be chosen based on the variance and mean of the likelihood scores from our posterior distribution parameter sets , using the ``leave-one-out" cross validation information criterion (LOO-IC). More variable and worse likelihood scores are penalized with this information criterion.

## Model list to start with
We can start with a quick list of the simplest 24 models that take combinatorix of overdispersion (none, sex-agnostic, sex-specific), excess zeroes (zero-hurdle, zero inflation), parental mortality (included, excluded), and sex-specific mean adjustment (included, excluded) into account. From that list of $3 \times 2 \times 2 \times 2 = 24$, we can more quickly see the worst models that are unreasonable to continue with and we can make further adjustments to the models considered realistic enough by LOO-IC (e.g. sex-specific excess zeroes, time-specific effects, interactions of the above terms).
