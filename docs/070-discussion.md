


# Discussion {#discussion}

## Model selection is not always the goal

Building a model motivated by a set of principles and domain expertise should be the preferred way of performing an analysis. The next important principle is model comparison, especially in terms of predictive inference. One model also doesn't always work for everything. In the course of building a model that is just complex enough to answer questions about age and temporal recalibration, I mentioned that intermediate models could be used to answer questions about average effects at different levels. For purely predictive inference, there is also the possibility of Bayesian model averaging (BMA) and other ensemble methods.

## Data Cleaning and Reproducibility

Data doesn't always come in a nice tidy^[Tidy data is described by Hadley Wickham] format, and I had the pleasure of turning the raw experimental data into a clean data set that is ready for modeling. Sometimes the process is quick and straight forward, but other times, like with this psychometric data, it takes more effort and clever techniques. There is academic value in describing the steps I took up front to reduce the headache later.

To begin, there is a strong push in recent years for reproducible data science. Scientific methods and results should be able to be replicated by other researchers, and part of that includes being able to replicate the process that takes the raw data and produces the tidy data that is ready for analysis. Tidy data is described by @wickham2014tidy and can be summed up by three principles

1. Each variable forms a column
2. Each observation forms a row
3. Each type of observational unit forms a table

One problem I have come across and have been guilty of in the past is having data in a spread sheet, modifying it, and then having no way of recovering the original data. Spread sheets are a convenient way to organize, transform, and lightly analyze date, but problems can quickly arise unless there is a good backup/snapshot system in place. Data is immutable^[Mutability in computer science is the property of a data structure where its contents can be modified in place. Immutability means that the object cannot be modified without first making a copy.], or at least that is the mindset that researchers must adopt in order to have truly reproducible workflows. The raw data that is collected or produced by a measurement device should never be modified without first being copied, even if for trivial reasons such as correcting a spelling mistake^[If a change is made to the raw data, it should be carefully documented and reversible].

To begin the data cleaning journey, I'll introduce the directory system that I had been given to work with. Each task is separated into its own folder, and within each folder is a subdirectory of age groups.

```
RecalibrationData
├── ParticipantAgeSex.xlsx
├── Audiovisual
│   ├── MiddleAge
│   ├── Older
│   └── Young
├── Duration
│   ├── MiddleAge
│   ├── Older
│   └── Young
├── Sensorimotor
│   ├── MiddleAge
│   ├── Older
│   └── Young
└── Visual
    ├── MiddleAge
    ├── Older
    └── Young
```

Within each age group subdirectory are the subdirectories for each subject named by their initials which then contain the experimental data in Matlab files.

```
├── MiddleAge
│   ├── BT
│   │   ├── BTadapt1__MAT.mat
│   │   ├── ...
│   ├── ...
├── Older
│   ├── BB
│   │   ├── BBadapt1__MAT.mat
│   │   ├── ...
│   ├── ...
└── Young
    ├── AC
    │   ├── ACadapt1__MAT.mat
    │   ├── ...
    ├── ...
```

At this point, the data appears manageable, there is information contained in the directory structure such as task, age group, and initials, and file name contains information about the experimental block. There is also an excel file that I was later given that contains more subject information like age and sex, though that information is not used in the model. The columns of the Matlab file depends on the task, but generally contains an SOA value and a response, but no column or row name information - that was provided by the researcher who collected the data. 

The next thing I did was to create a table of metadata - information extracted from the directory structure and file names combined with the the subject data and the file path. Regular expressions can be used to extract patterns from a string. With a list of all Matlab files within the `RecalibrationData` folder, I tried to extract the task, age group, initials, and block using the expression

```
"^(\\w+)/(\\w+)/(\\w+)/[A-Z]{2,3}_*[A-Z]*(adapt[0-9]|baseline[0-9]*).*"
```

Breaking it apart, the `^(\\w+)/` matches any word characters at the start and before the next slash. Since the directory structure is `Task/AgeGroup/Subject/file.mat`, the regular expression should match three words between slashes. The file name generally follows the pattern of `Initials__block#__MAT.mat`, so `[A-Z]{2,3}_*[A-Z]*` should match the initials, and `(adapt[0-9]|baseline[0-9]*)` should match the block (baseline or adapt). This method works for $536$ of the $580$ individual records. For the ones it failed, it was generally do to misspellings or irregular capitalizing of "baseline" and "adapt".


```r
table(feat_typ[,4])
#> 
#>     adapt1     adapt2     adapt3   baseline  baseline2 baseline43 
#>        165        162         43        162          3          1
```


```r
table(feat_atyp[,4])
#> 
#>                         adapat1      adapation1    Adaptation_1    Adaptation_2 
#>               1               1               1               9               7 
#>     Adaptation1     adaptation2     Adaptation2     adaptation3          adpat1 
#>               2               1               3               1               1 
#>          adpat2          adpat3        baseline        Baseline Visual_Baseline 
#>               4               1               1              10               1
```

Since there is only a handful of irregular block names, they can be dealt with a separate regular expression that properly extracts the block information. Other challenges in cleaning the data include the handling of subjects with the same initials. This becomes a problem because filtering by a subject's initials is not guaranteed to return a unique subject. Furthermore there are two middle age subjects with the same initials of "JM", so one was also identified with their sex "JM_F". The solution is to create a unique identifier (labeled as SID) that is a combination of age group, sex, and initials. For an experiment identifier (labeled as RID), the task and block were prepended to the SID. Each of these IDs uniquely identify the subjects and their experimental records making it easier to filter and search.



```r
glimpse(features)
#> Rows: 580
#> Columns: 8
#> $ rid       <fct> av-post1-M-f-CC, av-post1-M-f-DB, av-post1-M-f-HG, av-post1…
#> $ sid       <fct> M-f-CC, M-f-DB, M-f-HG, M-f-JM, M-f-MS, M-f-SJF, M-f-TS, M-…
#> $ path      <chr> "Audiovisual/MiddleAge/CC/CCadapt1__MAT.mat", "Audiovisual/…
#> $ task      <chr> "audiovisual", "audiovisual", "audiovisual", "audiovisual",…
#> $ trial     <fct> post1, post1, post1, post1, post1, post1, post1, post1, pos…
#> $ age_group <fct> middle_age, middle_age, middle_age, middle_age, middle_age,…
#> $ age       <dbl> 39, 44, 41, 48, 49, 43, 47, 49, 49, 44, 43, 44, 48, 48, 50,…
#> $ sex       <fct> F, F, F, F, F, F, F, F, F, M, M, M, M, M, M, F, F, F, F, F,…
```


Then with the table of clean metadata, the task is simply to loop through each row, read the Matlab file given by `path`, add the unique ID as a column, and then join the experimental data with the metadata to create a data set that is ready for model fitting and data exploration. The full code used to generate the clean data is not yet available online, but can be shared with the committee.

The benefit of writing a script to generate the data is that others can look over my code and verify that it is doing what I intended for it to do, and I can go back to any step within the process to make changes if the need comes up. Another tool that contributed to the reproducibility is the version control management software, Git. With Git I can take a snapshot of the changes I make, and revert if necessary. This thesis is also hosted on Github, and the entire history of development can be viewed there.

## Developing a model

[Chapter 3](#workflow) details the deeper considerations that went into building a model, but doesn't tell the full story of struggles and setbacks I faced. I find that I learn more from others when they share what didn't work along with the final path that did work. There is knowledge to be gained in failed experiments, because then there is one more way to not do something, just like a failing outcome reduces the variance of the Beta distribution.

I knew that I wanted to apply Bayesian modeling techniques to the data, because it was something knew that I was learning. I tried using a classical GLM to first get a baseline understanding of the data, but the fact that some estimates for certain subjects failed due to complete separation reinforced my enthusiasm to employ non-classical techniques. My first Bayesian model was derived from @lee2014bayesian which used nested loops to iterate over subjects and SOA values. I felt that the data was stored in a complicated way and made it difficult to comprehend and extend.

Next I moved on to using `arm::bayesglm` to remove convergence issues, but was met with other limitations such as linear parameterization and lack of hierarchical modeling. The book Statistical Rethinking [@mcelreath2020statistical] was my first introduction to Bayesian multilevel modeling. His `rethinking` package accompanies the book, and offers a compact yet expressive syntax for models that get translated into a Stan model. A model with age group and block can be written using `rethinking::ulam` as



```r
rethinking::ulam(alist(
  k ~ binomial_logit(n, p),
  p = exp(b + bG[G] + bT[trt]) * (x - (a + aG[G] + aT[trt])),
  a ~ normal(0, 0.06),
  aG[G] ~ normal(0, sd_aG),
  aT[trt] ~ normal(0, sd_aT),
  b ~ normal(3, 1),
  bG[G] ~ normal(0, sd_bG),
  bT[trt] ~ normal(0, sd_bT),
  c(sd_aG, sd_aT, sd_bG, sd_bT) ~ half_cauchy(0, 5)
), data = df, chains = 4, cores = 4, log_lik = TRUE)
```


During my time learning about multilevel models, I tried writing my own package that generates a Stan program based on R formula syntax. At the time I didn't fully understand the concepts of no-pooling, complete pooling, and partial pooling, and my package was plagued by the same lack of flexibility that `rstanarm` and `brms` have. In fact I learned that `brms` and `rstanarm` already did what I was trying to do after I had already started making my library, but it was a fun learning and programming experience. The fossilized remains of my attempt can be viewed on github.

I also tried using `lme4`, `rstanarm`, and `brms`, and learned more about the concepts of fixed and random effects. It was around this time that I noticed that parameterization can have a significant affect on the efficiency of a model and the inferential power of the estimated parameters. When fitting a classical model, there is little difference in estimating `a + bx` vs. `d(x - c)` since the latter can just be expanded as `-cd + dx` which is essentially the same as the first parameterization, but there is a practical difference in the interpretation of the parameters. The second parameterization implies that there is a dependence among the parameters that can be factored out. In the context of psychometric functions, there is a stronger connection between PSS and `c` and the JND and `d`. This parameterization made it easier to specify priors and also increased the model efficiency. Since only `rethinking` and `Stan` allow for arbitrary parameterization, I left the others behind.

I finally arrived at a model that worked well, but learned that using a binary indicator variable for the treatment comes with the assumption of higher uncertainty for one of the conditions. The linear model that I arrived at is displayed in equation \@ref(eq:badlinearmodel).


\begin{equation}
  \theta = \exp(\beta + \beta_G +(\beta_T + \beta_{TG})\times trt) \left[x - (\alpha + \alpha_G + (\alpha_T + \alpha_{TG})\times trt)\right]
  (\#eq:badlinearmodel)
\end{equation}


Using an indicator variable in this fashion also introduced an interaction effect into the model that I almost did not account for after I switched to using a factor variable. Interaction effects between factors is handled by creating a new factor that is essentially the cross-product of other factor variables. E.g. for factor variables $x$ and $y$

$$
x = \begin{bmatrix}
a \\
b \\
c
\end{bmatrix}, y =  \begin{bmatrix}
i \\
j
\end{bmatrix}\Longrightarrow x\times y = 
\begin{bmatrix}
ai & aj \\
bi & bj \\
ci & cj
\end{bmatrix}
$$

The final round of reparameterization came in the form of adopting non-centered parameterization for more efficient models. To us, $Z \sim N(0, 1^2);\quad X = 3 + 2Z$ is the same as $X \sim N(3, 2^2)$, but to a computer the process of sampling from $X$ can be more difficult than sampling from $Z$ (discussed in [chapter 4](#model-checking)).