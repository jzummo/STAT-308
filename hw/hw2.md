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

# Question 6

Show by simulation that the eigenvalues of the truncated covariance matrix of a stochastic process, $\Gamma_n$, densely approximate the range of the associated spectral density function as $n \rightarrow \infty$. Set $n = 1000$ and use an AR(1) process to verify.

### Solution

For our example, we simulate a Gaussian AR(1) procecss with $\phi = 0.5$. The commented lines prove the claim.

```{r}
n <- 1000
phi <- 0.5

acvs <- phi**(0:(n-1)) * (1/(1 - phi**2))
  
Gamma <- toeplitz(acvs)

eigenvalues <- eigen(Gamma)$values
print(range(eigenvalues))
# [1] 0.4444454 3.9999215

spectral_density <- function(omega) {
  1 / (1 - 2 * phi * cos(omega) + phi**2)
}

curve(spectral_density, from = 0, to = 2 * pi)
print(range(curve(spectral_density, from = 0, to = 2 * pi)$y))
# [1] 0.4444444 4.0000000
```



# Question 7

Let $\tilde \Gamma_n$ be the kernelized covariance matrix estimate with $j, h$ entry equal to
$$
\hat \gamma_{j-h} \one[|j - h| \leq k]
$$
Suppose we are working with an AR(1) process $X_i = \rho X_{i-1} + \epsilon_i$. Set $k = \rho \log n$. Prove that $\|\tilde\Gamma_n - \Gamma_n\| \cp 0$ in the spectral norm as $n \rightarrow \infty$.

### Solution

The AR(1) process has a representation as a linear process of the form
$$
X_i = \sum_{j=0}^\infty \rho^j \epsilon_{i-j}
$$
Thus the autocovariance is
$$
\gamma_k = \sum_{j=0}^\infty \rho^j \rho^{j+k} = \frac{\rho^{|k|}}{1 - \rho^2}.
$$
For fixed $\rho$, the important detail is that $\gamma_k \asymp \rho^{|k|}$. 

Next, we decompose the error into its respective estimation and approximation error:
$$
\tilde \Gamma_n - \Gamma_n = (\tilde \Gamma_n - \Gamma_n^{(k)}) + (\Gamma_n^{(k)} - \tilde \Gamma_n)
$$
where $\Gamma_n^{(k)}$ has $j, h$ entry equal to
$$
\gamma_{j-h}\one[|j - h| \leq k]
$$
By triangle inequality, it suffices to show that both errors converge to $0$ in the spectral norm.

The estimation error has $j, h$ component
$$
(\gamma_{j - h} - \hat \gamma_{j - h})\one[|j - h| \leq k]
$$
For $|x|_2 = 1$, we can bound
$$
x^\intercal(\tilde \Gamma_n - \Gamma_n^{(k)})x = \sum_{j,h}(\gamma_{j-h} - \hat \gamma_{j-h})\one[|j - h| \leq k] \\
= \sum_{|j - h| \leq k}(\gamma_{j-h} - \hat \gamma_{j-h}) \\
\leq (2k + 1)\max_{0 \leq \ell \leq k}|\gamma_\ell - \hat \gamma_\ell|
$$
Thus, it suffices to show that $\log n\cdot \max_{0 \leq \ell \leq k}|\gamma_\ell - \hat \gamma_\ell| \cp 0$.

We assume the $\epsilon_i$ are iid with mean $0$, variance $1$, and $E\epsilon_i^4 = \eta < \infty$, that $0 < \rho < 1$, and that $\hat\gamma_\ell = \frac{1}{n}\sum_{i=1}^{n-\ell} X_i X_{i+\ell}$.

Bias: Since $E\hat\gamma_\ell = \frac{n-\ell}{n}\gamma_\ell$,
$$
\log n \cdot \max_{0 \leq \ell \leq k}|E\hat\gamma_\ell - \gamma_\ell| = \log n \cdot \max_{0 \leq \ell \leq k}\frac{\ell}{n}|\gamma_\ell| \leq \frac{k \log n}{n(1-\rho^2)} = O\left(\frac{(\log n)^2}{n}\right) \rightarrow 0.
$$

Variance: Write $X_i = \sum_{a \geq 0}\psi_a \epsilon_{i-a}$ with $\psi_a = \rho^a$. The standard fourth-moment computation for linear processes gives, with $d = j - i$,
$$
\mathrm{Cov}(X_i X_{i+\ell}, X_j X_{j+\ell}) = \gamma_d^2 + \gamma_{d+\ell}\gamma_{d-\ell} + (\eta - 3)\sum_{m}\psi_m\psi_{m+\ell}\psi_{m+d}\psi_{m+d+\ell}.
$$
Using $|\gamma_r| \leq C\rho^{|r|}$ and $|d+\ell| + |d-\ell| \geq 2|d|$, each of the three terms is bounded by $C\rho^{2|d|}$, where $C$ does not depend on $\ell$. Hence, uniformly in $\ell$,
$$
\mathrm{Var}(\hat\gamma_\ell) \leq \frac{1}{n^2}\sum_{i,j=1}^{n} C\rho^{2|i-j|} \leq \frac{C'}{n}.
$$
By a union bound and Chebyshev's inequality, for any $\delta > 0$,
$$
P\left(\log n \cdot \max_{0 \leq \ell \leq k}|\hat\gamma_\ell - E\hat\gamma_\ell| > \delta\right) \leq \sum_{\ell=0}^{k}\frac{(\log n)^2\,\mathrm{Var}(\hat\gamma_\ell)}{\delta^2} \leq \frac{(k+1)(\log n)^2 C'}{n\delta^2} = O\left(\frac{(\log n)^3}{n}\right) \rightarrow 0.
$$
Combining the bias and variance bounds,
$$
\log n \cdot \max_{0 \leq \ell \leq k}|\hat\gamma_\ell - \gamma_\ell| \cp 0,
$$
and therefore $\|\tilde\Gamma_n - \Gamma_n^{(k)}\| \cp 0$.

Next, the approximation error has $j, h$ component
$$
\gamma_{j-h}\one[|j - h| > k]
$$
Let $x$ be a unit norm vector. Then
$$
x^\intercal (\Gamma_n^{(k)} - \tilde \Gamma_n)x = \sum_{|j - h| > k} \rho^{|j - h|}x_j x_k \\
\leq |\rho|^k \sum_{|j - h| > k} x_j x_k \\
\leq |\rho|^k |x|_2 \\
= |\rho|^k.
$$
Since $k \rightarrow \infty$ as $n \rightarrow \infty$, and $|\rho| < 1$, it follows that $|\rho|^k \rightarrow 0$ as $n \rightarrow \infty$, proving this approximation error goes to $0$ in the spectral norm.



# Question 8

Suppose $X_i \iid \N(0, 1)$. Prove that for any fixed positive integer $k$, the vector of scaled sample autocorrelations converges in distribution to a standard multivariate normal distribution as $n \to \infty$:
$$
\sqrt{n} \begin{pmatrix} \hat{\rho}_1 \\ \hat{\rho}_2 \\ \vdots \\ \hat{\rho}_k \end{pmatrix} \xrightarrow{d} \mathcal{N}(0, I_k),
$$

### Solution

Since the mean is known to be zero, the sample autocovariance and autocorrelation are
$$
\hat\gamma_h = \frac{1}{n}\sum_{i=1}^{n-h}X_i X_{i+h}, \qquad \hat\rho_h = \frac{\hat\gamma_h}{\hat\gamma_0}.
$$

We will assume the CLT for $m$-dependent stationary processes. We first prove a joint CLT for $(\hat\gamma_1, \dots, \hat\gamma_k)$ using the Cramér–Wold device, and then divide by $\hat\gamma_0$ and apply Slutsky's theorem.

By the Cramér–Wold device, it suffices to show that for every fixed $c = (c_1, \dots, c_k)^\intercal \in \R^k$,
$$
\sqrt{n}\sum_{h=1}^k c_h \hat\gamma_h \cd \N(0, c^\intercal c).
$$

The sums defining $\hat\gamma_h$ have different lengths, so we first replace them with sums over a common range. For each $h \leq k$,
$$
\sqrt{n}\,\hat\gamma_h - \frac{1}{\sqrt{n}}\sum_{i=1}^{n-k}X_i X_{i+h} = \frac{1}{\sqrt{n}}\sum_{i=n-k+1}^{n-h}X_i X_{i+h}.
$$

The right-hand side is a sum of at most $k$ terms, each $O_p(1)$, divided by $\sqrt{n}$, so it converges in probability to $0$. It therefore suffices to study
$$
\frac{1}{\sqrt{n}}\sum_{i=1}^{n-k} Y_i, \qquad Y_i = \sum_{h=1}^k c_h X_i X_{i+h}.
$$

The sequence $Y_i$ is stationary and $k$-dependent, since $Y_i$ depends only on $X_i, \dots, X_{i+k}$. It has mean $\E[Y_i] = \sum_h c_h \E[X_i]\E[X_{i+h}] = 0$.

We compute its long-run variance. For $h, h' \geq 1$,
$$
\E[X_i X_{i+h} X_j X_{j+h'}]
$$
is nonzero only if the four indices pair off into equal pairs. Since $i < i+h$ and $j < j+h'$, the only possible pairing is $i = j$ and $i + h = j + h'$, that is, $i = j$ and $h = h'$, in which case the expectation is $\E[X_i^2]\E[X_{i+h}^2] = 1$. The other pairing, $i = j + h'$ and $i + h = j$, would force $h = -h'$, which is impossible. Therefore
$$
\Cov[X_i X_{i+h}, X_j X_{j+h'}] = \one[i = j]\,\one[h = h'].
$$

It follows that
$$
\Var[Y_i] = \sum_{h,h'=1}^k c_h c_{h'}\one[h = h'] = \sum_{h=1}^k c_h^2 = c^\intercal c,
$$
and, for $m \geq 1$,
$$
\Cov[Y_i, Y_{i+m}] = \sum_{h,h'=1}^k c_h c_{h'}\one[i = i+m]\,\one[h = h'] = 0.
$$

So the long-run variance of $Y_i$ is $c^\intercal c + 2\sum_{m=1}^k 0 = c^\intercal c$. By the $m$-dependent CLT,
$$
\frac{1}{\sqrt{n}}\sum_{i=1}^{n-k} Y_i \cd \N(0, c^\intercal c).
$$

Since $c$ was arbitrary, the Cramér–Wold device gives
$$
\sqrt{n}\begin{pmatrix}\hat\gamma_1 \\ \vdots \\ \hat\gamma_k\end{pmatrix} \cd \N(0, I_k).
$$

It remains to normalize by $\hat\gamma_0$. By the weak law of large numbers, $\hat\gamma_0 = \frac{1}{n}\sum_{i=1}^n X_i^2 \cp \E[X_1^2] = 1$, so by the continuous mapping theorem $1/\hat\gamma_0 \cp 1$. Slutsky's theorem then gives
$$
\sqrt{n}\begin{pmatrix}\hat\rho_1 \\ \vdots \\ \hat\rho_k\end{pmatrix} = \frac{1}{\hat\gamma_0}\sqrt{n}\begin{pmatrix}\hat\gamma_1 \\ \vdots \\ \hat\gamma_k\end{pmatrix} \cd \N(0, I_k).
$$


# Question 9

This example is taken from a 1996 paper of Romano and Thombs. Let $X_i = \epsilon_i \epsilon_{i-1}$ where $\epsilon_i \iid \N(0, 1)$. Then $\E[X_i] = 0$ and
$$
\gamma_0 = \E[X_i^2] = \E[\epsilon_i^2]\E[\epsilon_{i-1}^2] = 1.
$$
And
$$
\gamma_1 = \E[\epsilon_i]\E[\epsilon_{i-1}^2]\E[\epsilon_{i-2}] = 0.
$$
Similarly, $\gamma_k = 0$ for $k \geq 1$. So $X_i$ is a white noise. Prove that
$$
\sqrt{n}(\hat \gamma_1 - \gamma_1) \Rightarrow \N(0, \sigma^2),
$$
where $\sigma^2 \neq 1$ (as we would expect under independence).

### Solution

We will assume the CLT for $m$-dependent stationary processes. Define $Y_i = X_i X_{i-1}$. Then $Y_i$ is $1$-dependent. We compute the long-run variance of $Y_i$ via
$$
\Var[Y_i] = \E[Y_i^2] \\
= \E[X_i^2 X_{i-1}^2] \\
= \E[\epsilon_i^2 \epsilon_{i-1}^2 \epsilon_{i-1}^2 \epsilon_{i-2}^2] \\
= \E[\epsilon_i^2]\E[\epsilon_{i-1}^4]\E[\epsilon_{i-2}^2] \\
= 3
$$

$$
\Cov[Y_i, Y_{i-1}] = \E[X_i X_{i-1}^2 X_{i-2}] \\
= \E[\epsilon_i \epsilon_{i-1} \epsilon_{i-1}^2\epsilon_{i-2}^2 \epsilon_{i-2} \epsilon_{i-3}] \\
= \E[\epsilon_i \epsilon_{i-1}^3\epsilon_{i-2}^3\epsilon_{i-3}] \\
= 0.
$$

So the long-run variance of $Y_i$ is $3 + 0 = 3$. By the $m$-dependent CLT,
$$
\frac{1}{\sqrt{n}}\sum_{i=1}^n Y_i \Rightarrow \N(0, 3)
$$
Since $Y_1 = o_p(\sqrt{n})$, we have the result:
$$
\frac{1}{\sqrt{n}}\sum_{i=2}^n Y_i \Rightarrow \N(0, 3).
$$


# Question 10

Consider the process $X_t = \rho |X_{t-1}| + \epsilon_t$ where $\epsilon_t \sim N(0, 1)$ and $|\rho| < 1$. Show that the pdf of the stationary marginal distribution is skew-normal with density proportional to
$$
\phi\left(x\sqrt{1 - \rho^2}\right)\Phi(\rho x)
$$

### Solution

Let $f(x)$ denote the relevant pdf. The distribution of $X_t$ conditional on $X_{t-1}$ is
$$
(X_t \mid X_{t-1}) \sim \N\left( \rho|X_{t-1}|, 1  \right).
$$
By the law of total probability, the stationary pdf of $X_t$ is
$$
f(x) = \int_\R f(y)\phi(x - \rho |y|)dy \\
= \int_0^\infty f(y) \phi(x - \rho y)dy + \int_{-\infty}^0 f(y) \phi(x + \rho y)dy \\
= \int_0^\infty f(y) \phi(x - \rho y)dy + \int_0^\infty f(-y) \phi(x - \rho y)dy \\
= \int_0^\infty \phi(x - \rho y)[f(y) + f(-y)]dy.
$$
It suffices to plug in our candidate pdf and verify the integral equation above. Observe that
$$
f(y) + f(-y) \propto \phi\left(y\sqrt{1 - \rho^2}\right)\Phi(\rho y) + \phi\left(-y\sqrt{1 - \rho^2}\right)\Phi(-\rho y) \\
= \phi\left(y\sqrt{1 - \rho^2}\right) \left( \Phi(\rho y) + \Phi(-\rho y) \right) \\
= \phi\left(y\sqrt{1 - \rho^2}\right).
$$
Thus, we must verify
$$
\phi\left(x\sqrt{1 - \rho^2}\right)\Phi(\rho x) \propto \int_0^\infty \phi(x - \rho y) \phi\left(y\sqrt{1 - \rho^2}\right) dy
$$
Dropping constants of proportionality,
$$
e^{-x^2(1 - \rho^2)/2}\Phi(\rho x) \propto \int_0^\infty e^{-(x - \rho y)^2/2}e^{-y^2(1 - \rho^2)}dy
$$

$$
\Phi(\rho x) \propto \int_0^\infty e^{x^2(1 - \rho^2)/2} e^{-(x - \rho y)^2/2}e^{-y^2(1 - \rho^2)/2}
$$

The exponential simplifies to
$$
\frac{-1}{2}\left( -(1 - \rho^2)x^2 + x^2 - 2 x \rho y + \rho^2 y^2 + y^2(1 - \rho^2) \right) \\
= \frac{-1}{2} \left( \rho^2 x^2 - 2x\rho y + y^2 \right) \\
= \frac{-1}{2} \left( \rho x - y \right)^2
$$
So the integral is
$$
\int_0^\infty \exp(-\left( \rho x - y \right)^2/2)dy
$$
Making the change of variables $\rho x - y = u$, the integral becomes
$$
-\int_{\rho x}^{-\infty} \exp(-u^2/2)du = \int_{-\infty}^{\rho x} \exp(-u^2/2)du = \Phi(\rho x).
$$

# Question 11

Extend the proof of GMC convergence to the case where $0 < p < 1$.

### Solution

The proof goes through as before when we substitute Minkowski's inequality for the elementary inequality $(a + b)^p \leq a^p + b^p$  as well as the completeness of $\ell_p$ as a metric space (not a normed space) when $0 < p < 1$.



# Question 12

Let $\rho$ denote the spectral radius and $\|\cdot\|_2$ the operator norm. 

- Design an example where the matrix $A$ satisfies $\rho(A) < 1$ but $|A|_2 > 1$. 
- Show that $\|A\|_2 < 1$ implies that $\rho(A) < 1$.

### Solution

For the first part, let
$$
A = 
\begin{bmatrix}
1 & 1 \\
0 & 1
\end{bmatrix}
$$
Since $A$ is upper triangular, the eigenvalues of $A$ are on the diagonal: $\{1, 1\}$. So $\rho(A) = 1$. However,
$$
A^\intercal A = 
\begin{bmatrix}
1 & 1 \\
1 & 2
\end{bmatrix}
$$
The characteristic equation is
$$
(\lambda - 2)(\lambda - 1) - 1 = 0
$$
The largest root is
$$
\lambda = \frac{3 + \sqrt{5}}{2}
$$
So the operator norm $\|A\|_2$ (equal to the largest singular value)
$$
\sqrt{\frac{3 + \sqrt{5}}{2}} > 1.
$$
For the second part, let $x$ be nonzero vector such that
$$
Ax = \rho(A) x.
$$
Taking the Euclidean norm on each side, we see
$$
\frac{|Ax|}{|x|} = \rho(A)
$$
By definition,
$$
\|A\|_2 = \sup_{x \neq 0} \frac{|Ax|}{|x|}
$$
It follows immediately that $\rho(A) \leq \|A\|_2$. 



