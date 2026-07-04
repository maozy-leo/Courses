# Inference on Linear Model



We first introduce a few denotations to clarify our problem setting.

Suppose $y=X\beta+\varepsilon$, where

- $y\in\R^n$ is the response vector
- $X$ is an $n\times p$ matrix of known constants(in this setting, we are referring to the fixed design, i.e. seeing $X$ as given constant matrix. Thus the following operator $\mathbb E$ or $Var$ are understood in the meaning of conditional expectation )
- $\beta\in\R^p$ is an unknown parameter vector
- $\varepsilon$ is a vector of random errors, satisfying $\mathbb E\varepsilon=0$

To do inference on linear model, we need to put some assumptions on the structure of the error $\varepsilon$.

## Gauss-Markov with Normal Assumption

A natural way is to assume that $\varepsilon$ follows the $i.i.d.$ Gaussian distribution $N(0,\sigma^2 I)$.

And in this setting we can directly derive many distributional results of the estimators, thus doing statistical inference and constructing confidence intervals.

### distribution of $\hat\sigma^2$ and $C\beta=d$

We first introduce the standard estimator of $\sigma^2$ and $C\beta$:
$$
\hat\sigma^2=\dfrac{Y'(I-P_X)Y}{n-\rank X}\\
C\hat\beta=C(X'X)^-X'Y
$$
The first result is $\dfrac{(n-\rank X)\hat\sigma^2}{\sigma^2}\sim\chi^2_{n-\rank X}$, which is given by the chi-square lemma, checking $(I-P_X)$ is idempotent.

Notice, the $\hat\sigma^2$ has a more known form in the sum of squared residuals:
$$
\begin{align}
(n-\rank X)\hat\sigma^2&=Y'(I-P_X)Y\\
&=((I-P_X)Y)'((I-P_X)Y) \\
&=(Y-X(X'X)^-XY)'(Y-X(X'X)^-XY)\\
&=(Y-X\hat\beta)'(Y-X\hat\beta)

\end{align}
$$
And $\hat\beta$ naturally follows the normal distribution $N(\beta,\sigma^2(X'X)^-)$ and $C\hat\beta \sim N(C\beta,\sigma^2C(X'X)^-C')$.

Another observation is that $\hat\sigma^2$ and $\hat\beta$ is independent, by checking $(X'X)^-X'(I-P_X)=0$.

With these distributional results, we can respectively construct the test statistics.

Suppose we want to test the null hypothesis $H_0:C\beta=d$. We assume that $C\in\R^{q\times p}$ is full rank in row, or we can simply get the equivalent null hypothesis with a simplified form.

Under the null hypothesis, $C\hat\beta-d\sim N(0,\sigma^2C(X'X)^-C')$.

Thus, we can construct a quadratic form
$$
(C\hat\beta-d)'[C(X'X)^-C']^{-1}(C\hat\beta-d)/\sigma^2\sim\chi^2_q
$$
Then, we get two independent quadratic form. Combining them gives the F-statistic we want:
$$
F=\dfrac{SSR/(n-\rank X)}{(C\hat\beta-d)'[C(X'X)^-C']^{-1}(C\hat\beta-d)/q}
$$
For a simplified version $c'\beta=d$, where $c\in\R^p$ is a vector instead of a matrix, we can construct t-statistic similarly.

The t-statistic is given by
$$
\begin{align}
t&=\dfrac{(c'\hat\beta-d)/(\sigma\sqrt{c'(X'X)^-c})}{\sqrt{\dfrac{\hat\sigma^2}{\sigma^2}}}\\
&=\dfrac{c'\hat\beta-d}{\hat\sigma\sqrt{c'(X'X)^-c}}

\end{align}
$$

## What if without normal assumption?

In this setting, the only information on structure of $\varepsilon$ we know is its homoskedestisty, i.e. $\mathbb E\varepsilon=0, Var(\varepsilon)=\sigma^2 I$.

Thus, to obtain more distributional information about the estimators, we need to consider the big sample properties.

### Asymptotic distribution of $\hat\beta$

We first have the decomposition from the definition of $\hat\beta$:
$$
\begin{align}
\hat\beta-\beta&=(X'X)^-X'Y-\beta\\
&=(X'X)^-X'X\beta-\beta+(X'X)^-X'\varepsilon\\
&=(X'X)^-X'\varepsilon
\end{align}
$$
We then have
$$
\begin{align}
\hat\beta-\beta&=(X'X/n)^-(X'\varepsilon/n)
\end{align}
$$
We consider the two terms respectively.

For the first term, we simply assume that the limit $\lim_{n\to\infty}\dfrac{X'X}{n}=Q$ exists and $Q$ is positive definite, to make sure that $(X'X/n)^{-1}$ is well-defined and non-zero.

For the second term, we want to use CLT to derive its limit distribution.
$$
\begin{align}
\dfrac{X'\varepsilon}{\sqrt n}&=\dfrac{1}{\sqrt n}\sum_{i=1}^nx_i\varepsilon_i
\end{align}
$$
To satisfy the Lindeberg condition, we put a few restrictions on the structure of $\varepsilon$ and $X$

- $\varepsilon_i$ are conditionally independent given $X$
- $\max_{1\le i\le n}\dfrac{||x_i||}{\sqrt n}\to0$

Use the multi-dimensional Lindeberg-Feller CLT, we have
$$
\dfrac{X'\varepsilon}{\sqrt n}\overset{d}{\to}N(0,\sigma^2 Q)
$$
Combining the estimation for the two terms gives the limit distribution of $\hat\beta$
$$
\sqrt n(\hat\beta-\beta)\overset d\to N(0,\sigma^2Q^{-1})
$$
A natural corollary is
$$
\sqrt n(C\hat\beta-C\beta)\overset d\to N(0,\sigma^2CQ^{-1}C')
$$

### Consistency of $\hat\sigma^2$

We still define $\hat\sigma^2$ as
$$
\hat\sigma^2=\dfrac{Y'(I-P_X)Y}{n-r}=\dfrac{\varepsilon'(I-P_X)\varepsilon}{n-r}
$$
In our assumption, $r=\rank X=p$ is constant, thus when $n\to\infty$, we only need to consider $\dfrac{\varepsilon'(I-P_X)\varepsilon}{n}=\dfrac1n\sum_{i=1}^n\varepsilon^2_i-\dfrac{1}{n}\varepsilon'P_X\varepsilon$.

The first term is in the form of sum of i.i.d. RV. Thus, with assumption $\mathbb E\varepsilon_i^4<\infty$, the SLLN gives $\dfrac{1}{n}\sum_{i=1}^n\varepsilon_i^2\overset{a.s.}\to\sigma^2$.

We now consider the second term.
$$
\begin{align}
\dfrac1n\varepsilon'P_X\varepsilon&=\dfrac{1}{n}\varepsilon'X(X'X)^{-1}X'\varepsilon\\
&=\left(\dfrac{X'\varepsilon}{n}\right)'\left( \dfrac{X'X}{n} \right)^{-1}\left(\dfrac{X'\varepsilon}{n}\right)

\end{align}
$$
From our previous discussion, we have
$$
\left(\dfrac{X'X}{n}\right)^{-1}\to Q\\
\left( \dfrac{X'\varepsilon}{n} \right)=O_p(\dfrac{1}{\sqrt n})=o_p(1)
$$
These results imply
$$
\dfrac1n\varepsilon'P_x\varepsilon=o_p(1)
$$
Thus, 
$$
\hat\sigma^2\overset p\to\sigma^2
$$

## Another view of the test $C\beta=d$

Let's begin with the blocked model, and for simplicity we still assume $\varepsilon\sim N(0,\sigma^2I)$.

Our model is written as
$$
Y=X_A\beta_A+X_B\beta_B+\varepsilon
$$
where
$$
X=\begin{pmatrix}
X_A,X_B
\end{pmatrix}
,\,\,\,\,
\beta=\begin{pmatrix}
\beta_A\\
\beta_B
\end{pmatrix}
$$
$X,\beta$ are partitioned so that the product $X_A\beta_A, X_B\beta_B$ are meaningful.

In this setting, we want to test the null hypothesis
$$
H_0:\beta_B=0
$$
Intuitively, $H_0$ means that covariates $X_B$ has no explanatory power of $Y$. Thus, the following two models are roughly the "same"
$$
Y=X_A\beta_A+\varepsilon_r\\
Y=X_A\beta_A+X_B\beta_B+\varepsilon_{ur}
$$
To be more accurate, they have almost the same explanatory power, i.e. the residual $\varepsilon_{ur}$ should not be much smaller than $\varepsilon_r$.

This enlightens us to construct the following F statistics
$$
F=\dfrac{(SSR_r-SSR_{ur})/(\rank(X)-\rank(X_A))}{SSR_{ur}/(n-\rank(X))}
$$
which can be interpreted as the ratio of compared with the restricted model, how much improvement can be gained if introducing new covariates $X_B$.

Clearly, under null hypothesis, there should be no improvement.
Now, we compute the exact distribution of the F-statistics.
We begin with the square sum of residuals of two models.
$$
\begin{align}
SSR_r&=(Y-X_A\hat\beta_A)'(Y-X_A\hat\beta_A)\\
&=Y'(I-P_{X_A})Y

\end{align}
$$

$$
\begin{align}
SSR_{ur}&=(Y-X\hat\beta)'(Y-X\hat\beta)\\
&=Y'(I-P_X)Y\\

\end{align}
$$

We now derive another form of $P_X$, so that we can get the reduced form of $P_X-P_{X_A}$.
To get a unique form, we assume that $X$ is full rank in column.
Consider the Schur complement of $(X'X)^{-1}$.

We can get its inverse in blocked form. (Detailed proof can be obtained in appendix)
$$
P_X=P_{X_A}+M_{X_A}X_B(X_B'M_{X_A}X_B)^{-1}X_B'M_{X_A}=P_{X_A}+P_{M_{X_A}X_B}
$$

<u>Remark</u>: This decomposition gives a direct sum decomposition
$$
Col(X)=Col(X_A)\oplus Col(M_{X_A}X_B)
$$

Thus,
$$
SSR_r-SSR_{ur}=Y'(P_X-P_{X_A})Y
$$
Thus,$\dfrac{SSR_r-SSR_{ur}}{\sigma^2}\sim\chi^2$ with degree of freedom $\rank(X)-\rank(X_A)$.
On the other hand, $\dfrac{SSR_{ur}}{\sigma^2}\sim\chi^2$ with degree of freedom $n-\rank(X)$.
$(I-P_X)(P_X-P_{X_A}=0)$ implies the independence of the two chi-square random variable.
Thus, $F\sim F(\rank(X)-\rank(X_A), n-\rank(X))$.



Back to our original problem, if we can convert our $H_0$ into $\tilde H_0: \tilde \beta=0$, we can use the above result to construct test statistics.
We begin with a simple identity
$$
Y=X\beta+\varepsilon=X(P_{C'}+I-P_{C'})\beta+\varepsilon=XC'(CC')^{-1}C\beta +X(I-P_{C'})\beta+\varepsilon
$$
To get a full rank design matrix, we rewrite $I-P_{C'}$ by its eigen vectors, which gives
$$
Y=XC'(CC')^{-1}\cdot C\beta+XU\cdot U'\beta+\varepsilon
$$
We then have
$$
\tilde Y=Y-XC'(CC')^{-1}d=XC'(CC')^{-1}\cdot (C\beta-d)+XU\cdot U'\beta+\varepsilon=Z_A\gamma_A+Z_B\gamma_B+\varepsilon
$$

which can be tested using our blockwise test method.

<u>Remark</u>: Both $XU$ or $X(I-P_{C'})$ are right choice for $Z_B$. The key point is that it should be full-ranked, representing the residual part after projection onto $Col(C')$.

## Generalization: what if without homoskedasticity?

In the previous part, we have introduced the statistics used for inference in the classical Gaussian-Markovian model.

In the following part, we will generalize our inference statistics to the hetoroskedasticity case, which is widely used in more circumstances.

We begin with the Aitken model
$$
Y=X\beta+\varepsilon
$$
where error $\varepsilon$ follows a general Gaussian distribution
$$
\varepsilon\sim N(0,\sigma^2V)
$$
where $\sigma^2$ is an unknown parameter and $V$ is a given positive definite matrix.

With a simple linear transformation $Z=V^{-1/2}Y$, the model is then in a form we are very familiar with
$$
Z=W\beta+e
$$
where $Z=V^{-1/2}Y, W=V^{-1/2}X, e=V^{-1/2}\varepsilon$. Notice, $e$ follows standard Gaussian distribution $e\sim N(0,\sigma^2)$.

Thus, all the inference results we get in the previous part, can be directly used in the Aitken model.

Specifically, for the null hypothesis
$$
H_0:C\beta=d
$$
We can easily derive the F-type statistics
$$
F=\dfrac{(C\hat\beta-d)'(C'(W'W)^-C)^{-1}(C\hat\beta-d)/\rank q}{Z'(I-P_W)Z/(n-\rank W)}
$$
where $W,Z$ are defined in the previous context, and $\hat\beta$ is given by $\hat\beta=(W'V^{-1}W)^{-1}W'V^{-1}Y$.



## Appendix

### Inverse of partitioned matrix
We want to find the inverse of the following blocked matrix.
To make notations clean, we obmit the assumption of invertibilities of certain matrices.
$$
\begin{pmatrix}
A&B\\C&D
\end{pmatrix}
$$
Perform elementary row&column transformations on the matrix to make it diagnol blockwise.
$$
\begin{pmatrix}
I&O\\
-CA^{-1}&I
\end{pmatrix}
\begin{pmatrix}
A&B\\C&D
\end{pmatrix}
=\begin{pmatrix}
A&B\\
O&D-CA^{-1}B
\end{pmatrix}
$$
Then
$$
\begin{pmatrix}
I&O\\
-CA^{-1}&I
\end{pmatrix}
\begin{pmatrix}
A&B\\C&D
\end{pmatrix}
\begin{pmatrix}
I&-A^{-1}B\\
O&I
\end{pmatrix}
=\begin{pmatrix}
A&O\\
O&D-CA^{-1}B
\end{pmatrix}
$$
Then we have
$$
\begin{pmatrix}
A&B\\C&D
\end{pmatrix}=

\begin{pmatrix}
I&O\\
CA^{-1}&I
\end{pmatrix}
\begin{pmatrix}
A&O\\
O&D-CA^{-1}B
\end{pmatrix}
\begin{pmatrix}
I&A^{-1}B\\
O&I
\end{pmatrix}
$$
Take inverse, we have
$$
\begin{pmatrix}
A&B\\C&D
\end{pmatrix}^{-1}=\begin{pmatrix}
I&-A^{-1}B\\
O&I
\end{pmatrix}
\begin{pmatrix}
A^{-1}&O\\
O&(D-CA^{-1}B)^{-1}
\end{pmatrix}
\begin{pmatrix}
I&O\\
-CA^{-1}&I
\end{pmatrix}
$$
Specifically, for $X=(X_1,X_2)$, inverse of $X'X$ is given by
$$

$$
