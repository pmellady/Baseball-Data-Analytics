Baseball Analytics Questions
================
Patrick Mellady

- [Problem 1](#problem-1)
  - [Model Set-up](#model-set-up)
  - [Data Considerations](#data-considerations)
  - [Implementation](#implementation)
  - [Results](#results)
- [Problem 2](#problem-2)
  - [Model Set-up](#model-set-up-1)
  - [Data Considerations](#data-considerations-1)
  - [Results](#results-1)
- [Appendix: Model Derivation for Question
  1](#appendix-model-derivation-for-question-1)
- [Appendix: Model Derivation for Question
  2](#appendix-model-derivation-for-question-2)
  - [The Multinomial is an Exponential Dispersion
    Family](#the-multinomial-is-an-exponential-dispersion-family)
  - [Determining the Mean Function in terms of
    $\theta$](#determining-the-mean-function-in-terms-of-theta)
  - [High Dimensional Form of the
    Model](#high-dimensional-form-of-the-model)
  - [MLE Estimation for $\beta$](#mle-estimation-for-beta)
  - [Adding a Ridge Penalty](#adding-a-ridge-penalty)

# Problem 1

*Develop a statistical model to predict the probability of a pitch being
called a strike, conditional on the batter not swinging.*

## Model Set-up

The data for this problem consist of $n=106,077$ observations of 23
variables. We have two binary variables `is_strike` and `is_swing`. We
will use these two binary variables to create a single $K=4$ dimension
multinomial vector, $Y_i$. We will use the following encoding:

- `is_strike`=1 and `is_swing`=0$\implies 1$
- `is_strike`=1 and `is_swing`=1$\implies 2$
- `is_strike`=0 and `is_swing`=1$\implies 3$
- `is_strike`=0 and `is_swing`=0$\implies 4$

With the above variable defined, we can proceed with a model definition.
To do this, we will introduce a latent Pólya-Gamma random variable for
each observation. This allows us to define the following hierarchical
model

$$
\begin{align*}
Y_i|B, b &\sim MN_4(1, \pi_i)\text{ where }\tilde\pi_i=f(\psi_i)\text{ and }\psi_i=X_iB+Z_ib\\
\omega_{ik}&\sim PG(n_{ik},0)\text{ for }i=1,2,3\cdots,n\text{ and }k=1,2,3,\\
B&\sim N(B_0, \Sigma_B)\\
b&\sim N(b_0, \Sigma_b)
\end{align*}
$$

Since our model is multinomial and we are working with a vectorized
version of the regression coefficients, as evidenced by the multivariate
normal prior on both $B$ and $b$, we must define $X_i$ as follows:

$$X_i=I_{K-1}\bigotimes x_i^T,\quad Z_i=I_{K-1}\bigotimes z_i^T$$

where $x_i$ and $z_i$ are the vector of fixed and random covariates for
observation i, respectively.

Additionally, the link function, $f$, is the stick breaking function.
This satisfies the following properties

$$
\begin{align*}
f(\psi_{i})=\frac{\exp(\psi_i)}{1+\exp(\psi_i)}=\tilde\pi_i\\
\tilde\pi_{ik}=\frac{\pi_{ik}}{1-\sum_{j<k}\pi_{ij}}
\end{align*}
$$

so that the multinomial probabilities are recovered from $\tilde\pi_i$
via $\pi_{ik}=\tilde\pi_{ik}\prod_{j<k}(1-\tilde\pi_{ij})$ for $k=1,2,3$
and $\pi_{i4}=\prod_{j=1}^3(1-\tilde\pi_{ij})$.

The above model yields the following conditional posterior distributions

$$
\begin{align*}
B|Y, b, \omega&\sim N((\sum_{i=1}^nX_i^T\Omega_iX_i+\Sigma_B^{-1})^{-1}(\sum_{i=1}^nX_i^T\Omega_i(\mu_i-Z_ib)+\Sigma_B^{-1}B_0), (\sum_{i=1}^nX_i^T\Omega_iX_i+\Sigma_B^{-1})^{-1})\\
b|Y, B, \omega&\sim N((\sum_{i=1}^nZ_i^T\Omega_iZ_i+\Sigma_b^{-1})^{-1}(\sum_{i=1}^nZ_i^T\Omega_i(\mu_i-X_iB)+\Sigma_b^{-1}b_0), (\sum_{i=1}^nZ_i^T\Omega_iZ_i+\Sigma_b^{-1})^{-1})\\
\omega_{ik}|Y, B, b &\sim PG(n_{ik}, \psi_{ik})
\end{align*}
$$

where

$$
\begin{align*}
\Omega_i=&\text{diag}(\omega_{ik}: k=1,2,3)\\
n_{ik}=&n_i-\sum_{j<k}Y_{ij}\\
\mu_{ik}=&\frac{1}{\omega_{ik}}(Y_{ik}-\frac{n_{ik}}{2})
\end{align*}
$$

Note that the definition of $\psi_i=X_iB+Z_ib$ allows for the use of
random effects in our model. Specifically, our data contains variables
indicating the pitcher, batter, catcher, and umpire. Setting the random
effect design matrix, $Z$, to the dummy encoded levels of the player
specific variables allows for each participant of a pitch event to
contribute a unique value to the probabilities of the outcomes.

Once we obtain posterior draws $B^{(s)}$ and $b^{(s)}$ for
$s=1,\dots,S$, we can calculate $\psi_i^{(s)}=X_iB^{(s)}+Z_ib^{(s)}$ for
each draw, transform it to a probability vector $\pi_i^{(s)}$ via $f$
and the stick-breaking construction, and average to obtain
$\hat\pi_i=\frac{1}{S}\sum_{s=1}^S\pi_i^{(s)}$. We may then calculate
the desired probability via:

$$
\hat P(\text{strike given no swing})=\frac{\hat\pi_{i1}}{\hat\pi_{i1}+\hat\pi_{i4}}
$$

## Data Considerations

The data for this problem present two glaring issues to be reconciled
before fitting: `spin_rate` has missing values and
`plate_location_x`/`plate_location_z` currently enter the design matrix
linearly.

To solve the missing data problem, we will impute the missing values of
`spin_rate` using random forests. This will allow for complex, nonlinear
relationships between the other variables and spin rate. This is done
through the `missForest` package.

To handle the nonlinearity of the strikezone, we will model plate
location using a tensor-product spline basis over horizontal and
vertical plate locations. This allows the model to learn a flexible
two-dimensional decision surface rather than imposing some sort of
smooth, parametric constraint.

## Implementation

In the appended code, helper functions were created to establish nicely
organized data and sample containers for use in the fitting process.
These functions are `build_matrices` and `initialize_samples`,
respectively.

The `build_matrices` function creates the design matrices for both the
fixed and random effects. It is important to note that, due to memory
constraints, these design matrices have **not** been operated on with
the Kronecker product yet. Lastly, this function returns two
$n\times(K-1)$ matrices corresponding to the $K-1$ dimensional
representation of the observations $Y_i$ and the values of $n_{ik}$.

The `initialize_samples` function creates a well organized list-of-lists
structure to hold the samples obtained at each iteration. For memory
constraints, we opt to omit the full trace of the latent $\omega$
variables, as they will be discarded after the samples are obtained
anyway.

The `MN_GIBB_2` function is what actually performs the Gibbs sampling
algorithm. This function takes the output of `build_matrices` as well as
the prior parameter specifications and returns the entire list-of-lists
of the generated samples.

The sampler was run for `nsamps`=5000 samples with the first 4000
samples being considered burn-in. The remaining 1000 samples were used
to generate estimates for $\psi_i$. The $\psi_i$ estimates were
transformed via the `stick_breaking_pi` function to obtain probability
vectors, which were then averaged and passed through the conditional
probability formula from the model set-up section to calculate the
desired conditional probability. This is all accomplished by the
`estimate_pi` function.

Lastly, the dataframe was modified to include the variable `is_event`.
This is a binary variable denoting whether or not a pitch was a strike
subsetted only to the observations where the batter did not swing. This
allows for fast computation of the model’s brier score, which is done
below.

## Results

We begin presentation of results by examining the convergence
diagnostics of the chain. Specifically, we observe the log likelihood
and Brier score trace plots across the run.

![](answers_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Next, we present the Brier score for our model along with Brier scores
for other comparative models. We compare our model to three others:

- A base model that predicts using the overall average strike
  probability
- A logistic regression model with fixed effects only
- A mixed effects logistic regression model

The models above are computed over the conditioning set where the batter
does not swing, while the multinomial model is fit over the entire data.
While this conditioning simplifies the problem and can lead to better
Brier scores, it uses only half of the data and does not account for
batter swing decisions. All Brier scores below are computed on the data
used to fit the models (in-sample).

| Model          | Brier Score |
|:---------------|------------:|
| Simple         |       0.215 |
| Logistic       |       0.056 |
| Logistic Mixed |       0.053 |
| Multinomial    |       0.063 |

The multinomial model produces a Brier score of 0.0630659, while the
simple mean model produced a Brier score of 0.2151697, and the logistic
model produced a Brier score of 0.0558078. This implies the multinomial
model outperforms the base model by 71% and underperforms the logistic
model by 13%.

Lastly, when comparing the multinomial model with a logistic mixed
effects model (fit using lme4), we find that the logistic GLMM produces
a Brier score of 0.05259 which makes the multinomial model 20% worse
than the logistic random effects model. While the logistic random
effects model produces a lower Brier score, the added ability to model
batter swing decisions presents a tradeoff. Depending on the context,
either model could be preferred. For example, the multinomial model
could be used to simulate all pitch-by-pitch outcomes for entire games,
in addition to providing high quality estimates for the estimated
conditional probability.

# Problem 2

*Build a projection for the strikeout rate for each pitcher in the
dataset for the 2025 season*

## Model Set-up

For this question, we will use a maximum likelihood approach. First, let
$Y_{ij}$ represent the observed response for player $i$ in year $j$ and
$S_{ij}$ represent the average stuff for player $i$ in year $j$. We can
then define the following model structure:

$$
Y_{ij}\sim MN_K(m_{ij}, \pi_{ij})\qquad\text{where }\pi_{ij}=\text{softmax}(\rho\cdot S_{i,j-1}+X_{ij}B)=\text{softmax}(\eta_{ij})
$$

where $\rho$ is a $(K-1)\times 1$ vector,
$X_{ij}=I_{K-1}\bigotimes x_{ij}^T$, and the linear predictor of the
$K^{th}$ category is fixed at 0 in the softmax. Under this framework, we
incorporate a pitcher-specific fixed effect via factor encoding in the
design matrix, $X_i$, and utilize stuff metrics via the one-year lag.
This specification allows us to incorporate the additional signal
captured in the stuff metrics without having to observe a players
current average stuff, which fits naturally with the goal of forecasting
rates for an unobserved year.

One issue to be found in this model is the sparse matrix layout of
$X_i$. With 200 pitchers across 5 years, the majority of entries of
$X_i$ will be zero, leading to estimation problems when computing the
MLE. To remedy this issue, we opt to implement a ridge regularization
when computing the coefficients $B$. To determine the optimal value of
the penalty parameter, we run a temporal cross-validation on our
dataset: we fit the model for the years 2020 through 2023, predict on
2024, and choose the penalty parameter that minimizes certain criteria.
In this case, the criteria for optimization is a combination of RMSE and
correlation between the observed and predicted rates.

## Data Considerations

As imported, the data currently shows simulated at-bats for 200 pitchers
across 5 years. This amounts to a dataframe with a very large number of
rows. By summarizing a pitcher’s entire season with aggregate counts for
strikeout and walks as well as mean stuff, we can collapse the data into
sufficient statistics for each pitcher-year. Additionally, this is a
complete dataset without any missingness, so imputation will not be a
consideration in our analysis.

Since this dataset is relatively sparse in covariates, we will make the
most use of the covariates at our disposal. Specifically, we will
implement higher-order parametric terms for both `year` and `age` to try
to capture the nonlinearities present in aging curves and overall league
dynamics.

## Results

Upon completing a cross-validation procedure, the model was fit to the
entire dataset. Producing the following metrics:

|           | Observed Rate | Predicted Rate |  RMSE |
|:----------|--------------:|---------------:|------:|
| Strikeout |         0.233 |          0.233 | 0.043 |
| Walk      |         0.082 |          0.081 | 0.028 |

![](answers_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Additionally, we can compare the above results with other models to
gauge the model’s fit. To do this, we fit two additional models:

- A purely autoregressive model of the form
  $Y_{ij}\sim MN_K(m_{ij}, \text{softmax}(A(Y_{i, j-1}/m_{i,j-1})))$
- A purely exogenous model of the form
  $Y_{ij}\sim MN_K(m_{ij}, \text{softmax}(X_{ij}B))$

The two models above yield the following results. For the exogenous only
model, we obtain:

|           | Observed Rate | Predicted Rate |  RMSE |
|:----------|--------------:|---------------:|------:|
| Strikeout |         0.233 |          0.233 | 0.043 |
| Walk      |         0.082 |          0.081 | 0.028 |

![](answers_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

and for the autoregressive only model, we get:

|           | Observed Rate | Predicted Rate |  RMSE |
|:----------|--------------:|---------------:|------:|
| Strikeout |         0.233 |          0.028 | 0.223 |
| Walk      |         0.082 |          0.002 | 0.092 |

![](answers_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

These results demonstrate that the exogenous covariates (age/pitcher)
are the most important factors to consider in projecting strikeout and
walk rates, as the purely exogenous model performs similarly to the
model with both exogenous and autoregressive covariates. Although the
performance is similar, the model with both variable types demonstrates
a slightly stronger correlation between the observed and predicted rates
when compared to the model with only exogenous covariates.

# Appendix: Model Derivation for Question 1

Recall that we assume the following

- $Y_i|B, b\sim MN_4(1, \pi_i)$ where $\tilde\pi_i=f(\psi_i)$ and
  $\psi_i=X_iB+Z_ib$
- $\omega_{ik}\sim PG(n_{ik},0)$ for $i=1,2,3\cdots,n$ and $k=1,2,3$
- $B\sim N(B_0, \Sigma_B)$
- $b\sim N(b_0, \Sigma_b)$

Now, see that, by breaking the multinomial mass function into a product
of independent binomial mass functions, we may write

$$
f(y_i|B, b)=\prod_{k=1}^{K-1}f(\psi_{ik})^{y_{ik}}(1-f(\psi_{ik}))^{n_{ik}-y_{ik}}=\prod_{k=1}^{K-1}\frac{\exp(\psi_{ik})^{y_{ik}}}{(1+\exp(\psi_{ik}))^{n_{ik}}}
$$

Using the fundamental Pólya-Gamma identity, we may write

$$
\frac{\exp(\psi_{ik})^{y_{ik}}}{(1+\exp(\psi_{ik}))^{n_{ik}}}=2^{-n_{ik}}\exp((y_{ik}-n_{ik}/2)\psi_{ik})\int_{0}^\infty \exp(-\omega_{ik}\psi_{ik}^2/2)p(\omega_{ik})d\omega_{ik}
$$

This implies that the full, joint model of the data and parameters is
given by

$$
P(Y, B, b, \omega)={\Big[}\prod_{i=1}^n\prod_{k=1}^{K-1}2^{-n_{ik}}\exp((y_{ik}-n_{ik}/2)\psi_{ik}-\omega_{ik}\psi_{ik}^2/2)p(\omega_{ik}){\Big]}p(B)p(b)
$$

We can complete the square in the exponential to obtain the following

$$
P(Y, B, b, \omega)={\Big[}\prod_{i=1}^n\prod_{k=1}^{K-1}2^{-n_{ik}}\exp(\frac{-1}{2/\omega_{ik}}(\psi_{ik}-\frac{Y_{ik}-n_{ik}/2}{\omega_{ik}})^2+\frac{(Y_{ik}-n_{ik}/2)^2}{2\omega_{ik}})p(\omega_{ik}){\Big]}p(B)p(b)
$$

By allowing $\Omega_i=\text{diag}(\omega_{ik}: k=1,2,\cdots,K-1)$ and
$\mu_{ik}=\frac{1}{\omega_{ik}}(Y_{ik}-\frac{n_{ik}}{2})$, we can use
vector algebra to rewrite this as

$$
P(Y, B, b, \omega)={\Big[}\prod_{i=1}^n\exp(\frac{-1}{2}(\psi_i-\mu_i)^T\Omega_i(\psi_i-\mu_i)+\frac{1}{2}\mu_i^T\Omega_i\mu_i)\prod_{k=1}^{K-1}2^{-n_{ik}}p(\omega_{ik}){\Big]}p(B)p(b)
$$

We can simplify further via

$$
P(Y, B, b, \omega)=\exp(\sum_{i=1}^n\frac{-1}{2}(\psi_i-\mu_i)^T\Omega_i(\psi_i-\mu_i)+\frac{1}{2}\mu_i^T\Omega_i\mu_i){\Big[}\prod_{i=1}^n\prod_{k=1}^{K-1}2^{-n_{ik}}p(\omega_{ik}){\Big]}p(B)p(b)
$$

In finding the conditional posterior for $B$, we have

$$
P(B|Y, b, \omega)=c\cdot P(Y, B, b, \omega)=c_1\cdot\exp(\sum_{i=1}^n\frac{-1}{2}(\psi_i-\mu_i)^T\Omega_i(\psi_i-\mu_i))p(B)
$$

And by Normal-Normal conjugacy, we obtain

$$
B|Y, b, \omega\sim N((\sum_{i=1}^nX_i^T\Omega_iX_i+\Sigma_B^{-1})^{-1}(\sum_{i=1}^nX_i^T\Omega_i(\mu_i-Z_ib)+\Sigma_B^{-1}B_0), (\sum_{i=1}^nX_i^T\Omega_iX_i+\Sigma_B^{-1})^{-1})
$$

Similarly, for $b$

$$
b|Y, B, \omega\sim N((\sum_{i=1}^nZ_i^T\Omega_iZ_i+\Sigma_b^{-1})^{-1}(\sum_{i=1}^nZ_i^T\Omega_i(\mu_i-X_iB)+\Sigma_b^{-1}b_0), (\sum_{i=1}^nZ_i^T\Omega_iZ_i+\Sigma_b^{-1})^{-1})
$$

Lastly, for $\omega_{ik}$, we have

$$
P(\omega|Y, B, b)=c\cdot P(Y, B, b, \omega)=c\cdot 2^{-n_{ik}}\exp(\frac{-1}{2/\omega_{ik}}(\psi_{ik}-\frac{Y_{ik}-n_{ik}/2}{\omega_{ik}})^2+\frac{(Y_{ik}-n_{ik}/2)^2}{2\omega_{ik}})p(\omega_{ik})
$$

which, by the exponential tilting property, gives us

$$
\omega_{ik}|Y, B, b\sim PG(n_{ik}, \psi_{ik})
$$

# Appendix: Model Derivation for Question 2

## The Multinomial is an Exponential Dispersion Family

Recall that if $Y\sim MN_k(m,\pi)$, then $Y$ has pmf given by

$$
f(y)=\frac{m!}{y_1!y_2!\cdots y_k!}\pi_1^{y_1}\pi_2^{y_2}\cdots\pi_k^{y_k}
$$

This can be written in the form of an exponential dispersion family. To
see this, we can exponentiate the natural log of the mass function to
obtain

$$
f(y)=\exp(\ln(f(y)))=\exp(\ln(m!)+\sum_{i=1}^ky_i\ln(\pi_i)-\ln(y_1!y_2!\cdots y_k!))
$$

We will convert the mass function to be parametrized in terms of the
first $k-1$ counts. To do this, note the $\pi_k=1-\sum_{i=1}^{k-1}\pi_i$
and $y_k=m-\sum_{i=1}^{k-1}y_i$. We can now write

$$
f(y)=\exp(\ln(m!)+\sum_{i=1}^{k-1}y_i\ln(\pi_i)+(m-\sum_{i=1}^{k-1}y_i)\ln(1-\sum_{i=1}^{k-1}\pi_i)-\ln(y_1!y_2!\cdots y_k!))
$$

Which can be simplified to

$$
f(y)=\exp(\sum_{i=1}^{k-1}y_i\ln(\frac{\pi_i}{1-\sum_{j=1}^{k-1}\pi_j})+m\ln(1-\sum_{j=1}^{k-1}\pi_j)+\ln(m!)-\ln(y_1!y_2!\cdots y_k!))
$$

If we let $Y=(y_1,y_2,\cdots,y_{k-1})$,
$\theta=(\ln(\frac{\pi_1}{1-\sum_{j=1}^{k-1}\pi_j}),\ln(\frac{\pi_2}{1-\sum_{j=1}^{k-1}\pi_j}),\cdots,\ln(\frac{\pi_{k-1}}{1-\sum_{j=1}^{k-1}\pi_j}))$,
$b(\theta)=-\ln(1-\sum_{j=1}^{k-1}\pi_j)$, and
$c(y,\phi)=\ln(m!)-\ln(y_1!y_2!\cdots y_k!)$, then we see that the
multinomial pmf has the form of an exponential dispersion family, so
standard frequentist GLM theory applies.

## Determining the Mean Function in terms of $\theta$

Consider that $\pi_k=1-\sum_{i=1}^{k-1}\pi_i$, so
$\ln(1-\sum_{i=1}^{k-1}\pi_i)=\ln(\pi_k)$. Now, we shall focus on a
single $\pi_i$ value. We asserted that
$\theta_i=\ln(\frac{\pi_i}{1-\sum_{j=1}^{k-1}\pi_j})=\ln(\frac{\pi_i}{\pi_k})$.
By exponentiating this equation, we have

$$
\pi_i=\pi_ke^{\theta_i}\implies 1=\sum_{i=1}^k\pi_i=\pi_k\sum_{i=1}^ke^{\theta_i}\implies\pi_k=\frac{1}{\sum_{i=1}^ke^{\theta_i}}
$$

hence

$$
b(\theta)=-\ln(1-\sum_{i=1}^{k-1}\pi_i)=-\ln(\pi_k)=\ln(\sum_{i=1}^ke^{\theta_i})
$$

This result further implies

$$
\pi_i=\frac{e^{\theta_i}}{\sum_{j=1}^ke^{\theta_j}}
$$

Further, see that for the $k^{th}$ category

$$
\theta_k=\ln(\frac{\pi_k}{\pi_k})=\ln(1)=0\implies e^{\theta_k}=1
$$

Which reduces our model to

$$
\pi_i=\frac{e^{\theta_i}}{1+\sum_{j=1}^{k-1}e^{\theta_j}},\qquad b(\theta)=\ln(1+\sum_{j=1}^{k-1}e^{\theta_j})
$$

## High Dimensional Form of the Model

In the usual GLM theory, our linear predictor is a scalar quantity:
$\eta_i=X_i^T\beta$. However, in the multinomial case, the linear
predictor is vector valued. Note that $g(\mu_i)=\eta_i$ and, since
$\mu_i$ is a $k-1$ dimensional vector and $g$ must be a one-to-one and
onto function, $\eta_i$ in this case must be a $k-1$ vector.

We can obtain this vector form of the linear predictor by expressing the
linear predictor as a matrix product

$$
\eta_i=x_i^T\begin{bmatrix} 
\beta_1 & \beta_2 & \cdots & \beta_{k-1}
\end{bmatrix}
$$

where $\beta_j$ is a $p_x\times1$ vector of coefficients corresponding
to each outcome category. In order to exploit the large sample theory
common to maximum likelihood estimation, must find another way to
express this equation, which we show below

$$
\eta_i=\begin{bmatrix}x_i^T & 0 & \cdots & 0\\
                      0 & x_i^T & \cdots & 0\\
                      \vdots & \vdots & \ddots & \vdots\\
                      0 & 0 & \cdots & x_i^T\end{bmatrix}\begin{bmatrix}\beta_1 \\ \beta_2 \\ \cdots \\ \beta_{k-1}\end{bmatrix}=X_iB
$$

where $X_i=I_{k-1}\bigotimes x_i^T$ and
$B=\begin{bmatrix}\beta_1^T&\beta_2^T&\cdots&\beta_{k-1}^T\end{bmatrix}^T$

## MLE Estimation for $\beta$

Using the above, we can write the log likelihood for a multinomial model
as

$$
l(B)=\sum_{i=1}^n{\Big[}\sum_{j=1}^{k-1}y_{ij}\theta_{ij}-m_i\ln(1+e^{\theta_{i1}}+\cdots+e^{\theta_{i,k-1}})+c(y_i,\phi){\Big]}
$$

Using vector notation, we can write this as

$$
l(B)=\sum_{i=1}^ny_i^T\theta_i-m_ib(\theta_i)+c(y_i,\phi)
$$

Now, see that we can write

$$
\nabla_B l(B)=\nabla_B\sum_{i=1}^ny_i^T\theta_i-m_ib(\theta_i)+c(y_i,\phi)=\sum_{i=1}^n\nabla_B(y_i^T\theta_i-m_ib(\theta_i))=\sum_{i=1}^n\nabla_Bl_i(B)
$$

where $l_i(B)=y_i^T\theta_i-m_ib(\theta_i)$. We can now use matrix
calculus to compute derivatives quickly for this model. First, using the
chain rule, we have

$$
\nabla_Bl_i(B)=\nabla_B(y_i^T\theta_i-m_ib(\theta_i))=\nabla_B\theta_i^T\nabla_{\theta_i}(y_i^T\theta_i-m_ib(\theta_i))
$$

Under the canonical link, $\theta_i=\eta_i$, so,
$\nabla_B\theta_i=\nabla_B\eta_i=X_i$. Further, we see that

$$
\nabla_{\theta_i}(y_i^T\theta_i-m_ib(\theta_i))=y_i-m_i\nabla_{\theta_i}b(\theta_i)
$$

Using the definition of $b(\theta_i)$ from the last section, we obtain

$$
\nabla_{\theta_i}b(\theta_i)=\begin{bmatrix}
\frac{e^{\theta_{i1}}}{1+\sum_{j=1}^{k-1}e^{\theta_{ij}}}\\
\frac{e^{\theta_{i2}}}{1+\sum_{j=1}^{k-1}e^{\theta_{ij}}}\\
\vdots\\
\frac{e^{\theta_{ik-1}}}{1+\sum_{j=1}^{k-1}e^{\theta_{ik-1}}}
\end{bmatrix}=\pi_i
$$

Thus, the gradient of the log likelihood is given by

$$
\nabla_Bl(B)=\sum_{i=1}^nX_i^T(y_i-m_i\pi_i)
$$

For the hessian matrix, we see that

$$
\nabla^2_Bl(B)=\nabla_B[\nabla_Bl(B)]=\nabla_B\sum_{i=1}^n\nabla_Bl_i(B)=\sum_{i=1}^n\nabla^2_Bl_i(B)
$$

we then can use the Bartlett identities to say

$$
E(\nabla^2_Bl_i(B))=-E(\nabla_Bl_i(B)[\nabla_Bl_i(B)]^T)
$$

This tells us that, via the Bartlett identities and the law of large
numbers, we have:

$$
E(\nabla^2_Bl(B))\approx-\sum_{i=1}^n\nabla_Bl_i(B)[\nabla_Bl_i(B)]^T=-\sum_{i=1}^nX_i^T(y_i-m_i\pi_i)(y_i-m_i\pi_i)^TX_i
$$

We can then employ the Fisher Scoring algorithm (Newton’s Method) to
find the MLE for our data.

## Adding a Ridge Penalty

To incorporate a ridge penalty in our model, we modify the objective we
are trying to optimize, specifically, we will determine the value of $B$
that maximizes the following quantity:

$$
l(B)-\frac{\lambda}{2}B^TB
$$

Since the derivative is a linear operator, the above calculations remain
relatively untouched, with the only changes we need to make being:

$$
\nabla_Bl(B)={\Big[}\sum_{i=1}^nX_i^T(y_i-m_i\pi_i){\Big]}-\lambda B
$$

and

$$
E(\nabla^2_Bl(B))\approx-\sum_{i=1}^n\nabla_Bl_i(B)[\nabla_Bl_i(B)]^T=-{\Big[}\sum_{i=1}^nX_i^T(y_i-m_i\pi_i)(y_i-m_i\pi_i)^TX_i{\Big]}-\lambda I
$$

which, again, we optimize the objective via Iteratively Reweighted Least
Squares (Newton’s Method).
