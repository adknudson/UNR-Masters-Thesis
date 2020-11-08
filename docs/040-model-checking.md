


# Model Fitting/Checking {#model-checking}

_Check your model before you wreck your model_

This chapter serves as the formal home of definitions and explanations of concepts relating to Markov Chain Monte Carlo (MCMC) and other diagnostic tools when working with Bayesian inference models. I touched on the physics of Hamiltonian Monte Carlo (HMC) and the diagnostic tools that come with it in the previous chapter, but now I will go into more detail.

## Fitting using HMC

Why do we need a sampler at all? Bayesian statistics and modeling stems from Bayes theorem (Equation \@ref(eq:bayesthm)). The prior $P(\theta)$ is some distribution over the parameter space and the likelihood $P(X | \theta)$ is the probability of an outcome in the sample space given a value in the parameter space. To keep things simple, we generally say that the posterior is proportional to the prior times the likelihood. Why proportional? The posterior distribution is a probability distribution, which means that the sum or integral over the parameter space must evaluate to one. Because of this constraint, the denominator in \@ref(eq:bayesthm) acts as a scale factor to ensure that the posterior is valid. 


\begin{equation}
  P(\theta | X) = \frac{P(X | \theta)\cdot P(\theta)}{\sum_i P(X | \theta_i)} =   \frac{P(X | \theta)\cdot P(\theta)}{\int_\Omega P(X | \theta)d\theta}
  (\#eq:bayesthm)
\end{equation}


For simple models, the posterior distribution can sometimes be evaluated analytically. An example of this is in _conjugate_ models, where the resulting posterior distribution is of the same type as the prior distribution, and an example of a conjugate model is the Beta distribution for inference about a proportion statistic. This is common in baseball for a player's batting average. I don't know a lot about baseball, but I know that hitting a baseball is a little less common than one in three swings, so _a priori_ I believe the probability of hitting a baseball is distributed as $\mathrm{Beta}(2, 5)$ because the expected value is $\approx 0.29$ and not a lot of weight is given to any particular value. Throughout a game I follow one player and he hits four balls and misses six - data that can be modeled as a Binomial observation. To figure out the posterior distribution for batting average, I use Bayes' theorem - _posterior is proportional to the prior times the likelihood_.


\begin{align*}
  P(\pi | y) &\propto P(y | \pi) \cdot P(\pi) \\
  &= {10 \choose 4}\pi^{4} (1-\pi)^{6} \cdot \frac{\Gamma(2+5)}{\Gamma(2)\Gamma(5)} \pi^{2-1}(1-\pi)^{5-1} \\
  &\propto \pi^{4+2-1}(1-\pi)^{6+5-1} \\
  &= \pi^{6-1}(1-\pi)^{11-1}
\end{align*}


The final line is the shape of a Beta distribution with parameters $6=2+4$ and $11=5+6$. The simple update rule is that for a prior $\mathrm{Beta}(a, b)$ and observed data with $y$ successes in $n$ observations, the posterior distribution is $\mathrm{Beta}(a + y, b + n - y)$. For the baseball player, the Bayesian estimate of his batting average is $6/(6+11) \approx 0.353$, but still with a good amount of uncertainty as shown in figure \@ref(fig:ch040-Teal-Monkey). 


<div class="figure" style="text-align: center">
<img src="040-model-checking_files/figure-html/ch040-Teal-Monkey-1.png" alt="After observing 4 hits in 10, the Beta(2,5) prior gets updated to become a Beta(6,11) posterior." width="85%" />
<p class="caption">(\#fig:ch040-Teal-Monkey)After observing 4 hits in 10, the Beta(2,5) prior gets updated to become a Beta(6,11) posterior.</p>
</div>


Conjugate models are great for simple observational data, but often it happens that the posterior distribution cannot be deduced from the model or that the integral in the denominator is complex or of a high dimension. In the former situation, the integral may not be possible to evaluate, and in the latter there may not be enough computational resources in the world to perform a simple grid approximation.

The solution is to use Markov Chain Monte Carlo (MCMC). The idea is that we can _draw samples_ from the posterior distribution in a way that samples proportional to the density. This sampling is a form of approximation to the area under the curve (i.e. an approximation to the denominator in \@ref(eq:bayesthm)). Rejection sampling [@gilks1992adaptive] and slice sampling [@neal2003slice] are basic methods for sampling from a target distribution, however they can often be inefficient^[Efficiency of a sampler is related to the proportion of proposal samples that get accepted.]. NUTS is a much more complex algorithm that can be compared to a physics simulation. A massless "particle" is flicked in a random direction with some amount of kinetic energy in a probability field, and is stopped randomly. The stopping point is the new proposal sample. The No U-Turn part means that when the algorithm detects that the particle is turning around, it will stop so as not to return to the starting position. This sampling scheme has a much higher rate of accepted samples, and also comes with many built-in diagnostic tools that let us know when the sampler is having trouble efficiently exploring the posterior. I'll talk more about these diagnostic tools throughout the remaining sections.

### Diagnostic Tools



#### Trace Plots

Trace plots are the first line of defense against misbehaved samplers. They are visual aids that let the practitioner asses the qualitative health of the chains, looking for properties such as autocorrelation, heteroskedacity, non-stationarity, and convergence. Healthy chains are _well-mixing_ and stationary. It's often better to run more chains during the model building process so that issues with mixing and convergence can be diagnosed sooner. Even one unhealthy chain can be indicative of a poorly specified model. The addition of more chains also contributes to the estimation of the Split $\hat{R}$ statistic, which I discuss in \@ref(split-r). Figure \@ref(fig:ch040-Brave-Moose) shows what a set of healthy chains looks like.

<div class="figure" style="text-align: center">
<img src="040-model-checking_files/figure-html/ch040-Brave-Moose-1.png" alt="An example of healthy chains." width="85%" />
<p class="caption">(\#fig:ch040-Brave-Moose)An example of healthy chains.</p>
</div>

There is a similar diagnostic plot called the rank histogram plot (or _trank_ plot for trace rank plot). @vehtari2020rank details the motivation for trank plots, but in short if the chains are all exploring the posterior efficiently, then the histograms will be similar and uniform. Figure \@ref(fig:ch040-Dog-Reborn) is from the same model as above but for the rank histogram.

<div class="figure" style="text-align: center">
<img src="040-model-checking_files/figure-html/ch040-Dog-Reborn-1.png" alt="A trank plot of healthy chains." width="85%" />
<p class="caption">(\#fig:ch040-Dog-Reborn)A trank plot of healthy chains.</p>
</div>

As the number of parameters in a model grows, it becomes exceedingly tedious to check the trace and trank plots of all parameters, and so numerical summaries are required to flag potential issues within the model.

#### $\hat{R}$ and Split $\hat{R}$ {#split-r}

The most common summary statistic for chain health is the potential scale reduction factor [@gelman1992inference] that measures the ratio of between chain variance and within chain variance. When the two have converged, the ratio is one. I've already shared examples of healthy chains which would also have healthy $\hat{R}$ values, but it's valuable to also share an example of a bad model. Below is the 8 Schools example [@gelman2013bayesian] which is a classical example for introducing Stan and testing the operating characteristics of a model. 



```r
schools_dat <- list(
  J = 8,
  y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
)
```





The initial starting parameters for this model are intentionally set to vary between $-10$ and $10$ (in contrast to the default range of $(-2, 2)$) and with only a few samples drawn in order to artificially drive up the split $\hat{R}$ statistic. The model is provided as supplementary code in the [appendix](#supplementary-code).



```r
fit_cp <- sampling(schools_mod_cp, data = schools_dat, refresh = 0,
                   iter = 50, init_r = 10, seed = 671254821)
```

Stan instantly warns about many different issues with this model, but the R-hat is the one of interest. The largest is $1.68$ which is incredibly large

<img src="040-model-checking_files/figure-html/ch040-Rocky-Test-1.png" width="85%" style="display: block; margin: auto;" />

These chains do not look good at all! Let's take a look at the $\hat{R}$ values and see if we can calculate one of the values manually.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-1)Split R-hat values from the 8 Schools example.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Parameter </th>
   <th style="text-align:right;"> Rhat </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mu </td>
   <td style="text-align:right;"> 1.234 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tau </td>
   <td style="text-align:right;"> 1.596 </td>
  </tr>
</tbody>
</table>

To calculate the (non split) $\hat{R}$, first calculate the between-chain variance, and then the average chain variance. For $M$ independent Markov chains, $\theta_m$, with $N$ samples each, the between-chain variance is

$$
B = \frac{N}{M-1}\sum_{m=1}^{M}\left(\bar{\theta}_m - \bar{\theta}\right)^2
$$

where

$$
\bar{\theta}_m = \frac{1}{N}\sum_{n=1}^{N}\theta_{m}^{(n)}
$$

and

$$
\bar{\theta} = \frac{1}{M}\sum_{m=1}^{M}\bar{\theta}_m
$$

The within-chain variance, $W$, is the variance averaged over all the chains.

$$
W = \frac{1}{M}\sum_{m=1}^{M} s_{m}^2
$$

where

$$
s_{m}^2 = \frac{1}{N-1}\sum_{n=1}^{N}\left(\theta_{m}^{(n)} - \bar{\theta}_m\right)^2
$$

The variance estimator is a weighted mixture of the within-chain and cross-chain variation

$$
\hat{var} = \frac{N-1}{N} W + \frac{1}{N} B
$$

and finally

$$
\hat{R} = \sqrt{\frac{\hat{var}}{W}}
$$

Here is the calculation in R


```r
param <- "mu"
theta <- p_cp[,,param]
N     <- nrow(theta)
M     <- ncol(theta)

theta_bar_m <- colMeans(theta)
theta_bar   <- mean(theta_bar_m)

B <- N / (M - 1) * sum((theta_bar_m - theta_bar)^2)
s_sq_m <- apply(theta, 2, var)

W <- mean(s_sq_m)
var_hat <- W * (N - 1) / N + B / N

(mu_Rhat <- sqrt(var_hat / W))
#> [1] 1.134
```

The $\hat{R}$ statistic is smaller than the split $\hat{R}$ value provided by Stan. This is a consequence of steadily increasing or decreasing chains. The split value does what it sounds like, and splits the chains in half and measures each half separately. In this way, the measure is more robust in detecting unhealthy chains. This also highlights the utility in using both visual and statistical tools to evaluate models.

#### Effective Sample Size

Samples from Markov Chains are typically autocorrelated, which can increase uncertainty of posterior estimates. I encountered this issue in the [second iteration](#iter2) of the model building process, and the solution I used was to reparameterize the model to avoid steep log-posterior densities - the benefit of reparameterization is conveyed by the ratio of effective sample size to actual sample size in figure \@ref(fig:ch040-Timely-Nitrogen). When the HMC algorithm is exploring difficult geometry, it can get stuck in regions of high densities, which means that there is more correlation between successive samples. 

<div class="figure" style="text-align: center">
<img src="040-model-checking_files/figure-html/ch040-Timely-Nitrogen-1.png" alt="Ratio of N_eff to actual sample size. Low ratios imply high autocorrelation which can be alleviated by reparameterizing the model or by thinning." width="85%" />
<p class="caption">(\#fig:ch040-Timely-Nitrogen)Ratio of N_eff to actual sample size. Low ratios imply high autocorrelation which can be alleviated by reparameterizing the model or by thinning.</p>
</div>

As the strength of autocorrelation generally decreases at larger lags, a simple prescription to decrease autocorrelation between samples and increase the effective sample size is to use _thinning_. Thinning means saving every $k^{th}$ sample and throwing the rest away. If one desired to have 2000 posterior draws, it could be done in two of many possible ways

- Generate 2000 draws after warmup and save all of them
- Generate 10,000 draws after warmup and save every $5^{th}$ sample. 

Both will produce 2000 samples, but the method using thinning will have less autocorrelation and a higher effective number of samples. Though it should be noted that generating 10,000 draws and saving all of them will have a higher number of effective samples than the second method with thinning, so thinning should only be favored to save memory.

#### Divergent Transitions {#divergent-transitions}

Unlike the previous tools for algorithmic faithfulness which can be used for any MCMC sampler, information about divergent transitions is intrinsic to Hamiltonian Monte Carlo. Recall that the HMC and NUTS algorithm can be imagined as a physics simulation of a particle in a potential energy field, and a random momentum is imparted on the particle. The sum of the potential energy and the kinetic energy of the system is called the Hamiltonian, and is conserved along the trajectory of the particle [@stanref]. The path that the particle takes is a discrete approximation to the actual path where the position of the particle is updated in small steps called _leapfrog steps_ (see @leimkuhler2004simulating for a detailed explanation of the leapfrog algorithm). A divergent transition happens when the simulated trajectory is far from the true trajectory as measured by the Hamiltonian.

A few divergent transitions is not indicative of a poorly performing model, and often divergent transitions can be reduced by reducing the step size and increasing the adapt delta parameter. On the other hand, a bad model may never be improved just by tweaking some parameters. This is the folk theorem of statistical computing - if there is a problem with the sampling, blame the model, not the algorithm.

Divergent transitions are never saved in the posterior samples, but they are saved internally to the Stan fit object and can be compared against good samples. Sometimes this can give insight into which parameters and which regions of the posterior the divergent transitions are coming from.

<div class="figure" style="text-align: center">
<img src="040-model-checking_files/figure-html/ch040-Hot-Locomotive-1.png" alt="Divergent transitions highlighted for some parameters from the second iteration model. Divergent transitions tend to occur when both the hierarchical variance terms are near zero." width="85%" />
<p class="caption">(\#fig:ch040-Hot-Locomotive)Divergent transitions highlighted for some parameters from the second iteration model. Divergent transitions tend to occur when both the hierarchical variance terms are near zero.</p>
</div>

## Prior Predictive Checks

I used prior predictive checks in the first iteration of the model to establish a few things pertaining to model adequacy and computational faithfulness. The first reason is to ensure that the selected priors do not put too much mass in completely implausible regions (such as really large JND estimates). Data simulated from the priors can also be used to check that the software works. When you have the exact priors that were used to generate the data, the fitting algorithm should be able to accurately recover the priors.



- transition to posterior predictive checks chapter
- fig 10 in for posterior predictive @gabry2019visualization