# Predictive Inference {#predictive-inferences}

_All models are wrong but some are useful_

The above quote is from George Box, and it is a popular quote that statisticians like to throw around^[I am one of them]. All models are wrong because it is nearly impossible to account for the minutiae of every process that contributes to an observed phenomenon, and often trying to results in poorer performing models. 

*why is predictive performance the right model selection/comparison criteria*

- idea of "geocentric" models: wrong models that still predict well  
- notions overfitting/underfitting:
- more parameters leads to better in-sample fit  
- a prefect fit to data is always possible  
- but predicts poorly (overfit)  
- underfitting fails to capture the *regular* features of the data (why regularizing priors are important)  


I think you covered this already in Ch. 1 and 2 but here is more thoughts:
The PI's predictive philosophy has evolved to prefer this reference model approach.
Early on statisticians are usually taught to prefer *parsimony* or simple models.
The idea is that this guards against *overfitting* and also boosts power to detect *statistically significant* effects.

Also computation limitations made small models preferable.
But in modern statistical learning, we tend to include all relevant data with elaborate probabilitistc structures.

The idea is to include all the data with the aim of squeezing all predictive ability from the data points.

- not sure where this goes, but make sure you say that 1 model is not sufficient, we need a collection (or series/sequence) of models. that is why we need to fit models fast in `stan`/HMC

transitional sentence: given that we want to compare models (and possibly select), how to quantifying


*Quantifying predictive performance*

- log posterior predictive (more below) and information theory (if you want to talk about that at all)
- cross-validation, loo, WAIC
- and estimates of loo. loo psis 
- @vehtari2017practical

*some notes from my grant posterior*. rewrite this for your glm based model.
Given a model $M$ with posterior predictive distribution $p( \tilde{T} | \tilde{x}, D$ for a new survival time $\tilde{T}$ with observed data $D$ with feature vector $\tilde{x}$. 
We evaluate predictive performance using the **logarithm of the predictive density (LPD)** evaluated pointwise at the actual observation $( \tilde{t}, \tilde{x}, M)$ [@Peltola2014; @Piironen2017b].
LPD is a proper scoring rule and measures both the **calibration** and **sharpness** of the predictive distribution [@Gneiting2007].
With omit technical definitions of these concepts, but loosely calibration means the statistical consistency between the predictive distribution and the observations (errors on the order).
Sharpness, on the other hand, refers to how concentrated the predictive posterior (how precisely forecasted).
Typically we don't have the analytic form of the predictive posterior, so instead we use $J$ MCMC draws to approximate the LPD [@Peltola2014]:

\begin{equation}
	LPD(M) \approx \frac{1}{J} \Sigma_{j=1}^{J} log p( \tilde{t} | \tilde{x}, D, \theta^{(j)} ),
\end{equation}

where $\theta^{(j)}$ is the posterior parameter vector from the $j$th posterior sample. 

Further we'll like a metric of general predictive performance and so compute the average over $n$ data points:

<!--
\begin{equation}
	MLPD(M) = \frac{1}{n} \Sigma_{i=1}^{n} log p( \tilde{t} | \tilde{x_i}, D, M ),
\end{equation}
-->

Further, we'd like to compare the MLPD value of a model $M$ and another model $M^*$ (possibly a reference model or competing model):

<!--
\begin{equation}
	\Delta MLPD(M) = MLPD(M) - MLPD(M^*)
\end{equation}
-->

A negative difference in $\Delta MLPD$ for Model $M$ compared to a reference Model ($M^*$) means worse performance for the model while a positive difference indicates better prediction.
We assess the uncertainty in the difference using Bayesian bootstrap [@Rubin1981] samples of $\Delta MLPD$ between model $M$ and $M^*$: