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

# Question 13

Let $\epsilon_i$ be a sequence of iid Rademacher random variables and let $X_i = \sum_{j=1}^\infty a_j \epsilon_{i - j}$ be a linear process.

- Assuming $a_j = j^{-1}$, prove that the partial sum process $S_n$ satisfies the following CLT: $\frac{S_n}{\sqrt{n} \log n} \Rightarrow \N(0, \sigma^2)$.
- Assuming $a_j = j^{-\beta}$ for $1/2 < \beta < 1$, prove that $\frac{S_n}{n^{3/2 - \beta}} \Rightarrow \N(0, \sigma^2)$. 

### Solution

$$
S_n = \sum_{i=1}^n \sum_{j=1}^\infty a_j \epsilon_{i - j} \\
= \sum_{t=1}^{n-1} \left(\sum_{i=1}^{n-t}a_i\right) \epsilon_t + \sum_{t=-\infty}^0 \left(\sum_{i=1}^n a_{i-t}\right)\epsilon_t
$$

Thus,
$$
\Var[S_n] = \sum_{t=1}^{n-1} \left(\sum_{i=1}^{n-t}a_i\right)^2 + \sum_{t=-\infty}^0 \left(\sum_{i=1}^n a_{i-t}\right)^2 \\
= \sum_{t=1}^{n-1} \left(\sum_{i=1}^{n-t}a_i\right)^2 + \sum_{t=-\infty}^0 \left(\sum_{i=1-t}^{n-t} a_i\right)^2
$$
For our first choice of coefficients,
$$
\Var[S_n] = \sum_{t=1}^{n-1} \left(\sum_{i=1}^{n-t}i^{-1} \right)^2 + \sum_{t=-\infty}^0 \left(\sum_{i=1-t}^{n-t} i^{-1}\right)^2
$$
By the integral comparison test,
$$
\sum_{i=1}^{n-t} i^{-1} \asymp \int_1^{n-t} \frac{dx}{x} = \log(n - t)
$$

$$
\sum_{i=1-t}^{n-t} i^{-1} \asymp \log(n - t) - \log(1 - t)
$$

Thus,
$$
\Var[S_n] \asymp \sum_{t=1}^{n-1} \log(n - t)^2 + \sum_{t=-\infty}^0 (\log(n - t) - \log(1 - t))^2
$$
We need to show that this is on the order of $n \log^2(n)$ for the CLT.

# TODO



# Question 14

Let $X_i = \sum_{j=0}^\infty a_j \varepsilon_{i-j}$ be a linear process with coefficients decaying as $a_i \sim i^{-\beta}$, where $\{\varepsilon_t\}$ are i.i.d. innovations with mean zero and finite fourth moments. Let CLT 1 denote the CLT for the sample mean $S_n/\sqrt{n}$, and CLT 2 for the sample variance. Using the predictive dependence measure $\theta_i$, prove that CLT 1 requires $\beta > 1$, while CLT 2 requires $\beta > 3/4$.

### Solution

Assume WLOG that $\epsilon_t$ is standard Gaussian white noise. The definition of predictive dependence is
$$
\theta_i = \|P_0 X_i\|_2 = \|\E[X_i \mid \mathcal F_0] - \E[X_i \mid \mathcal F_{-1}]\|_2
$$
We have
$$
\E[X_i \mid \mathcal F_0] = \sum_{j = i}^\infty a_j \epsilon_{i - j} \qquad \E[X_i \mid \mathcal F_{-1}] = \sum_{j = i + 1}^\infty a_j \epsilon_{i - j}
$$
The difference is $a_i \epsilon_0$. So the predictive dependence measure is $\theta_i = \|a_i \epsilon_0\|_2 = |a_i|$, and summability requires
$$
\sum_{i=1}^\infty |a_i| < \infty \iff \sum_{i=1}^\infty i^{-\beta} < \infty \iff \beta > 1.
$$
The sample variance CLT stipulates that
$$
\frac{1}{\sqrt{n}}\sum_{i=1}^n X_i^2 \Rightarrow \N(0, \sigma^2)
$$
We need to analyze the predictive dependence measure of $X_i^2$.
$$
\E[X_i^2 \mid \mathcal F_0] = \E\left[ \left( \sum_{j=0}^\infty a_j \epsilon_{i - j} \right)^2 \Bigg\vert \mathcal{F}_0 \right] \\
= \sum_{j = 0}^{i - 1} a_j^2 + \left(\sum_{j=i}^\infty a_j \epsilon_{i - j} \right)^2
$$
This is because of conditional orthogonality:
$$
\E\left[ \sum_{j=0}^{i - 1} a_j \epsilon_{i-j} \sum_{j=i}^\infty a_j \epsilon_{i - j} \Bigg\vert \mathcal F_0 \right] = \E\left[ \sum_{j=0}^{i - 1} a_j \epsilon_{i-j} \right] \sum_{j=i}^\infty a_j \epsilon_{i - j} = 0.
$$
Likewise,
$$
\E[X_i^2 \mid \mathcal F_{-1}] = \sum_{j=0}^i a_j^2 + \left(\sum_{j=i+1}^\infty a_j \epsilon_{i - j} \right)^2
$$
The difference is
$$
\theta_i = \left\lVert a_i^2 +  \left(\sum_{j=i+1}^\infty a_j \epsilon_{i - j} \right)^2 - \left(\sum_{j=i}^\infty a_j \epsilon_{i - j} \right)^2  \right\rVert_2
$$
Using the sum of squares identity $a^2 - b^2 = (a - b)(a+b)$:
$$
\left(\sum_{j=i+1}^\infty a_j \epsilon_{i - j} \right)^2 - \left(\sum_{j=i}^\infty a_j \epsilon_{i - j} \right)^2 = a_i \epsilon_0\left(a_i \epsilon_0 + 2\sum_{j = i + 1}^\infty a_j \epsilon_{i - j} \right)
$$
So we need to compute
$$
\theta_i^2 = \Var\left[ a_i \epsilon_0\left(a_i \epsilon_0 + 2\sum_{j = i + 1}^\infty a_j \epsilon_{i - j} \right) \right] \\
= \Var\left[ a_i^2 \epsilon_0^2 + 2 a_i \epsilon_0 \sum_{j = i + 1}^\infty a_j \epsilon_{i - j} \right]
$$
We use the law of total variance conditioning on $\epsilon_0$. The conditional variance is
$$
\Var\left[ a_i^2 \epsilon_0^2 + 2a_i \epsilon_0 \sum_{j = i + 1}^\infty a_j \epsilon_{i - j} \bigg\vert \epsilon_0 \right] = 4a_i^2 \epsilon_0^2 \Var\left[\sum_{j = i + 1}^\infty a_j \epsilon_{i - j}\right] \\
= 4a_i^2 \epsilon_0^2 \sum_{j = i + 1}^\infty a_j^2
$$
The expected conditional variance is thus $4a_i^2 \sum_{j = i + 1}^\infty a_j^2$. The conditional expectation is
$$
\E\left[ a_i^2 \epsilon_0^2 + 2a_i \epsilon_0 \sum_{j = i + 1}^\infty a_j \epsilon_{i - j} \bigg\vert \epsilon_0 \right] = a_i^2 \epsilon_0^2
$$


The variance of the conditional expectation is thus $2a_i^4$. Combining,
$$
\theta_i^2 = 2a_i^4 + 4a_i^2 \sum_{j = i + 1}^\infty a_j^2 \\
\asymp 2i^{-4\beta} + 4i^{-2\beta}\sum_{j = i + 1}^\infty j^{-2\beta}
$$
To find the asymptotic order of the tail sum, we use an integral comparison:
$$
\sum_{j=i+1}^\infty j^{-2\beta} \asymp \int_i^\infty x^{-2\beta}dx = \frac{i^{1-2\beta}}{2\beta - 1}
$$
Thus,
$$
\theta_i^2 \asymp i^{-4\beta} + i^{-2\beta} i^{1 - 2\beta} = i^{-4\beta} + i^{1 - 4\beta}
$$
So $\theta_i \asymp i^{1/2 - 2\beta}$. Summable physical dependence requires
$$
\sum_{i=1}^\infty \theta_i < \infty \iff \sum_{i=1}^\infty i^{1/2 - 2\beta} < \infty \iff 1/2 - 2\beta < -1
$$
Rearranging this inequality gives $3/4 < \beta$, the necessary condition for the sample variance process to satisfy a CLT.