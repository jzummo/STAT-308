$$
\newcommand{\E}{\mathbb{E}}
\newcommand{\Pr}{\mathbb{P}}
\newcommand{\Var}{\operatorname{Var}}
\newcommand{\Cov}{\operatorname{Cov}}
\newcommand{\Corr}{\operatorname{Corr}}
\newcommand{\R}{\mathbb{R}}

\newcommand{\ed}{\overset{d}{=}}
\newcommand{\cd}{\overset{d}{\rightarrow}}
\newcommand{\cp}{\overset{P}{\rightarrow}}
\newcommand{\as}{\text{a.s.}}
\newcommand{\cas}{\overset{\as}{\rightarrow}}
\newcommand{\cellone}{\overset{L^1}{\rightarrow}}

\newcommand{\st}{\text{ s.t. }}
\newcommand{\rv}{\text{ r.v. }}
\newcommand{\io}{\text{ i.o. }}
\newcommand{\iid}{\overset{\mathrm{iid}}{\sim}}

\newcommand{\one}{\mathbf{1}}

\newcommand{\argmax}{\operatorname*{arg\,max}}
\newcommand{\argmin}{\operatorname*{arg\,min}}

\newcommand{\iprod}[2]{\left\langle #1,#2\right\rangle}

\newcommand{\X}{\mathcal{X}}
\newcommand{\Y}{\mathcal{Y}}
\newcommand{\ind}{\perp\!\!\!\perp}

\newcommand{\N}{\mathcal{N}}
\newcommand{\Pois}{\operatorname{Poisson}}
\newcommand{\Exp}{\operatorname{Exponential}}
\newcommand{\Ber}{\operatorname{Bernoulli}}
\newcommand{\Bin}{\operatorname{Binomial}}
\newcommand{\Unif}{\operatorname{Uniform}}
\newcommand{\BetaDist}{\operatorname{Beta}}
\newcommand{\Gam}{\operatorname{Gamma}}

\newcommand{\lnorm}[2]{\left\|#1\right\|_{#2}}
\newcommand{\indc}[1]{\mathbf{1}_{\left\{#1\right\}}}
\newcommand{\norm}[1]{\left\|#1\right\|}
\newcommand{\hf}{1/2}

\newcommand{\wh}{\widehat}
\newcommand{\wt}{\widetilde}

\newcommand{\Norm}[1]{\left\|#1\right\|}
\newcommand{\Fnorm}[1]{\left\|#1\right\|_{\mathrm F}}
\newcommand{\fnorm}[1]{\left\|#1\right\|_{\mathrm F}}
\newcommand{\opnorm}[1]{\left\|#1\right\|_{\mathrm{op}}}
\newcommand{\nnorm}[1]{\left\|#1\right\|_{*}}

\newcommand{\Span}{\operatorname{span}}
\newcommand{\rank}{\operatorname{rank}}
\newcommand{\Tr}{\operatorname{Tr}}
\newcommand{\diag}{\operatorname{diag}}
\newcommand{\supp}{\operatorname{supp}}

\newcommand{\pth}[1]{\left(#1\right)}
\newcommand{\qth}[1]{\left[#1\right]}
\newcommand{\sth}[1]{\left\{#1\right\}}

\newcommand{\calL}{\mathcal{L}}
\newcommand{\TV}{\mathrm{TV}}
$$

# Question 1

Let $\epsilon_t$ be an iid sequence with $\E[\log^+ |\epsilon_t|] < \infty$. Prove that $X_t = \sum_{j=0}^\infty \alpha^j \epsilon_{t-j}$ converges almost surely for $|\alpha| < 1$ using Kolmogorov's three-series theorem. Here $\log^+(x) = \max(\log(x), 0)$ is the nonnegative logarithm.

### Solution

This theorem says that $\sum_{j=0}^\infty Y_j$ converges almost surely if, for any $c > 0$, letting $Y_j^c = Y_j \one[Y_j \leq c]$ be the truncation:
$$
\sum_{j=0}^\infty \Pr(|Y_j| > c)
$$

$$
\sum_{j=0}^\infty \E[Y_j^c] < \infty
$$

$$
\sum_{j=0}^\infty \Var[Y_j^c] < \infty
$$

We can take $c = 1$ for simplicity and $t = 0$. Set $Y_j = \alpha^j \epsilon_{-j}$.

##### (a)

$$
\sum_{j=0}^\infty \Pr(|\alpha^j||\epsilon_{-j}| > 1) \\
= \sum_{j=0}^\infty \Pr(|\epsilon_{-j}| > |\alpha|^{-j}) \\
\leq \sum_{j=0}^\infty \Pr( \log^+|\epsilon_{-j}| > \log (|\alpha|^{-j})) \\
= \sum_{j=0}^\infty \Pr( \log^+|\epsilon_0| > j\log |\alpha|^{-1}) \\
$$

Since $\E[\log^+|\epsilon_0|] < \infty$ and $\log^+ |\epsilon_0|$ is nonnegative, the layer cake representation gives us
$$
\E[\log^+ |\epsilon_0|] = \int_0^\infty \Pr(\log^+ |\epsilon_0| > t)dt \\
\propto \int_0^\infty \Pr(\log^+ |\epsilon_0| > t \log |\alpha|^{-1})dt,
$$
being finite by assumption. This integrand is monotone decreasing, so
$$
\sum_{j=1}^\infty \Pr(\log^+ |\epsilon_0| > j \log |\alpha|^{-1}) \leq \int_0^\infty \Pr(\log^+ |\epsilon_0| > t \log |\alpha|^{-1})dt < \infty
$$
Thus
$$
\sum_{j=0}^\infty \Pr(\log^+ |\epsilon_0| > j \log |\alpha|^{-1}) = \Pr(\log^+|\epsilon_0| > 0) + \sum_{j=1}^\infty \Pr(\log^+ |\epsilon_0| > j \log |\alpha|^{-1}) < \infty.
$$

##### (b)

We need to show that
$$
\E[Y_j^c] = \E[\alpha^j \epsilon_{-j} \one[\alpha^j \epsilon_{-j} \leq 1]]
$$

is summable over the nonnegative integers. Observe that
$$
\E[\alpha^j \epsilon_{-j} \one[\alpha^j \epsilon_{-j} \leq 1]] \leq \E[\one[\alpha^j \epsilon_{-j} \leq 1]] \\
= \Pr(\alpha^j \epsilon_{-j} \leq 1)
$$
Because of the rapid geometric decay,
$$
\Pr( \alpha^j \epsilon_{-j} \leq 1 \ \io ) = 0
$$
And these events $E_j = \{\alpha^j \epsilon_{-j} \leq 1\}$ are independent, so we can apply the second Borel-Cantelli lemma to deduce
$$
\sum_{j=0}^\infty \Pr( \alpha^j \epsilon_{-j} \leq 1 ) < \infty
$$

##### (c)

Similarly,
$$
\Var(Y_j^c) \leq \one[\alpha^{2j}\epsilon_{-j}^2 \leq 1]
$$
The same Borel-Cantelli argument thus proves that $\sum_{j=0}^\infty \Var[Y_j^c] < \infty$.



# Question 2

Let $X_1, \dots, X_n$ be drawn from a stationary stochastic process. Find the BLUE for estimating the common mean $\mu$ of the marginal distribution.

### Solution

Let $\Gamma_n$ denote the sample autocovariance matrix (using the $1/n$ normalization convention). The BLUE estimator can be written as
$$
\tilde \mu = \sum_{i=1}^n w_i X_i = w^\intercal X
$$
Unbiasedness requires that $\sum_{i=1}^n w_i = 1$. The variance of $\tilde \mu$ is $w^\intercal \Gamma_n w$. The Lagrangian for this linearly constrained quadratic program is thus
$$
w^\intercal \Gamma_n w - \lambda (w^\intercal \mathbf{1} - 1)
$$


 The stationarity conditions are
$$
2\Gamma_n w - \lambda \mathbf{1} = 0 \iff w = \frac{\lambda}{2}\Gamma_n^{-1}\mathbf{1}
$$
Since $\Gamma_n$ is positive semidefinite, this is a standard convex optimization problem with convex constraints, and duality applies. With the normalization constraint, the solution becomes
$$
w = \frac{\Gamma_n^{-1} \mathbf{1}}{\mathbf{1}^\intercal \Gamma_n^{-1} \mathbf{1}}
$$
Thus the BLUE is
$$
\tilde \mu = \frac{\mathbf{1}^\intercal \Gamma_n^{-1} X}{\mathbf{1}^\intercal \Gamma_n^{-1} \mathbf{1}}
$$

# Question 3

What is the variance of the BLUE estimator for a stationary process?

### Solution

With $w = \frac{\Gamma_n^{-1} \mathbf{1}}{\mathbf{1}^\intercal \Gamma_n^{-1} \mathbf{1}}$ as above, the variance is given by
$$
\frac{\one^\intercal \Gamma_n^{-1} \Gamma_n \Gamma_n^{-1}\one}{(\mathbf{1}^\intercal \Gamma_n^{-1} \mathbf{1})^2} \\
= \frac{1}{\one^\intercal \Gamma_n^{-1} \one}.
$$

# Question 4

Construct an example to prove that $\Gamma_n$ need not be consistent for $\Gamma$, in the sense that $\rho(\hat\Gamma_n - \Gamma_n) \not\rightarrow 0$ need not hold in general.

### Solution

We'll take the simplest case of Gaussian white noise for our example and $n = 300$ across $10$ simulations. The difference between the sample and true autocovariance matrix fluctuates like a random walk and does not tamper down to $0$ as expected from consistency.

```{r}
folder_path <- "~/Downloads"

for (b in 1:10) {
  file_name <- file.path(folder_path, paste0("my_plot_", b, ".png"))
  png(filename = file_name, width = 6, height = 4, units = "in", res = 300)
  
  results = c()
  x <- c()
  for (n in 1:300) {
    x <- append(x, rnorm(1))
    
    hatGamma <- acf(x, lag.max = n - 1, type = "covariance", plot = F)$acf[, 1, 1]
    hatGamma <- toeplitz(hatGamma)
    
    Gamma = diag(n)
    
    results <- append(results, norm(hatGamma - Gamma, type = "2"))
  }
  plot(results)
  
  dev.off()
}
```

![my_plot_10](/Users/jakezummo/Documents/STAT 308/images/my_plot_1.png)

![my_plot_2](/Users/jakezummo/Documents/STAT 308/images/my_plot_2.png)
![my_plot_3](/Users/jakezummo/Documents/STAT 308/images/my_plot_3.png)
![my_plot_4](/Users/jakezummo/Documents/STAT 308/images/my_plot_4.png)
![my_plot_5](/Users/jakezummo/Documents/STAT 308/images/my_plot_5.png)
![my_plot_6](/Users/jakezummo/Documents/STAT 308/images/my_plot_6.png)
![my_plot_7](/Users/jakezummo/Documents/STAT 308/images/my_plot_7.png)
![my_plot_8](/Users/jakezummo/Documents/STAT 308/images/my_plot_8.png)
![my_plot_9](/Users/jakezummo/Documents/STAT 308/images/my_plot_9.png)
![my_plot_10](/Users/jakezummo/Documents/STAT 308/images/my_plot_10.png)



# Question 5

Let $(X_0, X_k)^\top$ be a zero-mean bivariate Gaussian random vector defined by $(X_0, X_k)^\top \sim \N(0, \Sigma)$. Prove that
$$
\Corr(X_0, X_k)^2 = \Corr(X_0^2, X_k^2).
$$


 ### Solution

We can assume without loss of generality that $\Var(X_0) = \Var(X_k) = 1$ because of the identity $\Corr(|c|X, \cdot) = \Corr(X, \cdot)$. Let $\rho = \Corr(X_0, X_k)$. So we have
$$
\Sigma = 
\begin{bmatrix}
1 & \rho \\
\rho & 1
\end{bmatrix}
$$
We can orthognalize this as
$$
\Sigma = 
\left(\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 & 1 \\
1 & -1
\end{bmatrix}
\right)
\begin{bmatrix}
1 + \rho & 0 \\
0 & 1 - \rho\\
\end{bmatrix}
\left(\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 & 1 \\
1 & -1
\end{bmatrix}
\right)
$$
So
$$
\Sigma^{-1/2} = 
\left(\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 & 1 \\
1 & -1
\end{bmatrix}
\right)
\begin{bmatrix}
(1 + \rho)^{-1/2} & 0 \\
0 & (1 - \rho)^{-1/2}\\
\end{bmatrix}
\left(\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 & 1 \\
1 & -1
\end{bmatrix}
\right)
$$
The point is that
$$
X_0 - X_k \ind X_0 + X_k
$$
(which we could have derived any number of ways). So
$$
(X_0 - X_k)^2 \ind (X_0 + X_k)^2
$$

$$
\E[(X_0 - X_k)^2 (X_0 + X_k)^2]  = \E[(X_0 - X_k)^2] \E[(X_0 + X_k)^2]
$$

The RHS is
$$
(\E[X_0^2] - 2\E[X_0 X_k] + \E[X_k^2])(\E[X_0^2] + 2\E[X_0 X_k] + \E[X_k^2]) \\
= (2 - 2\rho)(2 + 2\rho) \\
= 4 - 4\rho^2
$$
And the LHS is
$$
\E[(X_0^2 - X_k^2)^2] \\
= \E[X_0^4 - 2X_0^2X_k^2 + X_k^4] \\
= 6 - 2\E[X_0^2X_k^2]
$$
Equating them yields
$$
\E[X_0^2 X_k^2] = (3 - 2) + 2\rho = 1 + 2\rho^2
$$
Now,
$$
\Corr(X_0^2, X_k^2) \\
= \frac{\Cov(X_0^2, X_k^2)}{2} \\
= \frac{\E[(X_0^2 - 1)(X_k^2 - 1)]}{2} \\
= \frac{\E[X_0^2 X_k^2] - \E[X_0^2] - \E[X_k^2] + 1}{2} \\
= \frac{1 + 2\rho^2 - 1}{2} \\
= \rho^2,
$$
which is what we wanted to prove.

