# Problem 1

*Develop a statistical model to predict the probability of a pitch being called a strike, conditional on the batter not swinging.*

## Model Set-up
The data for this problem consist of $n=106,077$ observations of 23 variables. We have two binary variables `is_strike` and `is_swing`. We will use these two binary variables to create a single $K=4$ dimension multinomial vector, $Y_i$. We will use the following encoding:


* `is_strike`=1 and `is_swing`=0$\implies 1$
* `is_strike`=1 and `is_swing`=1$\implies 2$
* `is_strike`=0 and `is_swing`=1$\implies 3$
* `is_strike`=0 and `is_swing`=0$\implies 4$

With the above variable defined, we can proceed with a model definition. To do this, we will introduce a latent Pólya-Gamma random variable for each observation. This allows us to define the following hierarchical model

$$
\begin{align*}
Y_i|B, b&\sim MN_4(1, \pi_i)\text{ where }\tilde\pi_i=f(\psi_i)\text{ and }\psi_i=X_iB+Z_ib\\
\omega_{ik}&\sim PG(1,0)\text{ for }i=1,2,3\cdots,n\text{ and }k=1,2,3,\\
B&\sim N(B_0, \Sigma_B)\\
b&\sim N(b_0, \Sigma_b)
\end{align*}
$$

Since our model is multinomial and we are working with a vectorized version of the regression coefficients, as evidenced by the multivariate normal prior on both $B$ and $b$, we must define $X_i$ as follows:

$$
X_i=I_{K-1}\bigotimes x_i^T,\quad Z_i=I_{K-1}\bigotimes z_i^T
$$

where $x_i$ and $z_i$ are the vector of fixed and random covariates for observation i, respectively.

Additionally, the link function, $f$, is the stick breaking function. This satisfies the following properties

$$
\begin{align*}
f(\psi_{i})=\frac{\exp(\psi_i)}{1+\exp(\psi_i)}=\tilde\pi_i\\
\tilde\pi_{ik}=\frac{\pi_{ik}}{1-\sum_{j<k}\pi_{ij}}
\end{align*}
$$

so that the multinomial probabilities are recovered from $\tilde\pi_i$ via $\pi_{ik}=\tilde\pi_{ik}\prod_{j<k}(1-\tilde\pi_{ij})$ for $k=1,2,3$ and $\pi_{i4}=\prod_{j=1}^3(1-\tilde\pi_{ij})$.

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

Note that the definition of $\psi_i=X_iB+Z_ib$ allows for the use of random effects in our model. Specifically, our data contains variables indicating the pitcher, batter, catcher, and umpire. Setting the random effect design matrix, $Z$, to the dummy encoded levels of the player specific variables allows for each participant of a pitch event to contribute a unique value to the probabilities of the outcomes.

Once we obtain posterior draws $B^{(s)}$ and $b^{(s)}$ for $s=1,\dots,S$, we can calculate $\psi_i^{(s)}=X_iB^{(s)}+Z_ib^{(s)}$ for each draw, transform it to a probability vector $\pi_i^{(s)}$ via $f$ and the stick-breaking construction, and average to obtain $\hat\pi_i=\frac{1}{S}\sum_{s=1}^S\pi_i^{(s)}$. We may then calculate the desired probability via:

$$
\hat P(\text{strike given no swing})=\frac{\hat\pi_{i1}}{\hat\pi_{i1}+\hat\pi_{i4}}
$$

## Data Considerations
The data for this problem present two glaring issues to be reconciled before fitting: `spin_rate` has missing values and `plate_location_x`/`plate_location_z` currently enter the design matrix linearly.

To solve the missing data problem, we will impute the missing values of `spin_rate` using random forests. This will allow for complex, nonlinear relationships between the other variables and spin rate. This is done through the `missForest` package.

To handle the nonlinearity of the strikezone, we will model plate location using a tensor-product spline basis over horizontal and vertical plate locations. This allows the model to learn a flexible two-dimensional decision surface rather than imposing some sort of smooth, parametric constraint.

## Implementation
In the appended code, helper functions were created to establish nicely organized data and sample containers for use in the fitting process. These functions are `build_matrices` and `initialize_samples`, respectively.

The `build_matrices` function creates the design matrices for both the fixed and random effects. It is important to note that, due to memory constraints, these design matrices have **not** been operated on with the Kronecker product yet. Lastly, this function returns two $n\times(K-1)$ matrices corresponding to the $K-1$ dimensional representation of the observations $Y_i$ and the values of $n_{ik}$.

The `initialize_samples` function creates a well organized list-of-lists structure to hold the samples obtained at each iteration. For memory constraints, we opt to omit the full trace of the latent $\omega$ variables, as they will be discarded after the samples are obtained anyway.

The `MN_GIBB_2` function is what actually performs the Gibbs sampling algorithm. This function takes the output of `build_matrices` as well as the prior parameter specifications and returns the entire list-of-lists of the generated samples.

The sampler was run for `nsamps`=5000 samples with the first 4000 samples being considered burn-in. The remaining 1000 samples were used to generate estimates for $\psi_i$. The $\psi_i$ estimates were transformed via the `stick_breaking_pi` function to obtain probability vectors, which were then averaged and passed through the conditional probability formula from the model set-up section to calculate the desired conditional probability. This is all accomplished by the `estimate_pi` function.

Lastly, the dataframe was modified to include the variable `is_event`. This is a binary variable denoting whether or not a pitch was a strike subsetted only to the observations where the batter did not swing. This allows for fast computation of the model's brier score, which is done below.

## Results
We begin presentation of results by examining the convergence diagnostics of the chain. Specifically, we observe the log likelihood and Brier score trace plots across the run.

```{r, echo=FALSE}
res<-readRDS("MN_GIBBS_samples.rds")

iters<-1:length(res$log_lik)

p1<-ggplot()+geom_line(aes(x=iters, y=res$log_lik))+
  labs(title="Log Likelihood Trace Plot")+
  xlab("Iteration Number")+ylab("Log Likelihood")
p2<-ggplot()+geom_line(aes(x=iters, y=res$bs))+
  labs(title="Brier Score Trace Plot")+
  xlab("Iteration Number")+ylab("Brier Score")

p1 | p2
```

Next, we present the Brier score for our model along with Brier scores for other comparative models. We compare our model to three others:

* A base model that predicts using the overall average strike probability
* A logistic regression model with fixed effects only
* A mixed effects logistic regression model

The models above are computed over the conditioning set where the batter does not swing, while the multinomial model is fit over the entire data. While this conditioning simplifies the problem and can lead to better Brier scores, it uses only half of the data and does not account for batter swing decisions. All Brier scores below are computed on the same data used to fit the models (in-sample).

```{r, echo=FALSE}
pitch_data<-read.csv("pitch_data_output.csv")

## Calculate the Brier score
brier<-mean((pitch_data$p_hat - pitch_data$is_event)^2, na.rm=TRUE)

brier_simple<-mean(
  (mean(pitch_data$is_event, na.rm=TRUE)-pitch_data$is_event)^2, na.rm=TRUE
)

X<-pitch_data[,c("is_event", "plate_location_x", "plate_location_z", "rel_speed", "spin_rate", 
            "induced_vert_break", "horizontal_break")]
X<-X %>% filter(!is.na(is_event))
X$plate_location_x_2<-X$plate_location_x^2
X$plate_location_z_2<-X$plate_location_z^2
X$plate_location_xz<-X$plate_location_x*X$plate_location_z

mod<-glm(is_event~., data=X, family=binomial)
brier_logistic<-mean((predict(mod, newdata = X, type="response")-X$is_event)^2, na.rm=TRUE)

no_swing<-read.csv("for_glmer.csv")
fit<-readRDS("glmer_fit.rds")

brier_glmm<-mean((predict(fit, newdata=no_swing, type="response")-no_swing$is_event)^2)

result_df<-data.frame(Model=c("Simple", "Logistic", "Logistic Mixed", "Multinomial"),
                      `Brier Score`=round(c(brier_simple, brier_logistic, brier_glmm, brier),3),
                      check.names=FALSE)

kable(result_df)
```

The multinomial model produces a Brier score of `r brier`, while the simple mean model produced a Brier score of `r brier_simple`,  and the logistic model produced a Brier score of `r brier_logistic`. This implies the multinomial model outperforms the base model by `r paste0(100*(1-round(brier/brier_simple, 2)),"%")` and the logistic model by `r paste0(100*(1-round(brier/brier_logistic, 2)),"%")`.

Lastly, when comparing the multinomial model with a logistic mixed effects model (fit using lme4), we find that the logistic GLMM produces a Brier score of `r brier_glmm` which makes the multinomial model `r paste0(100*(round(brier/brier_glmm, 2)-1),"%")` worse than the logistic random effects model. While the logistic random effects model produces a lower Brier score, the added ability to model batter swing decisions presents a tradeoff. Depending on the context, either model could be preferred. For example, the multinomial model could be used to simulate all pitch-by-pitch outcomes for entire games, in addition to providing high quality estimates for the estimated conditional probability.
