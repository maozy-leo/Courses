# Maximum Likelihood Estimation

Suppose the observed sample is

$$
X_1,\dots,X_n
$$

with joint density or probability mass function

$$
p(x_1,\dots,x_n\mid \theta),
$$

where $\theta\in\Theta\subseteq\mathbb R^k$ is unknown. The basic idea of maximum likelihood estimation is simple: among all distributions in the parametric family, choose the parameter value under which the observed sample is most likely to occur.

The likelihood function is

$$
L_n(\theta)=p(x_1,\dots,x_n\mid \theta).
$$

The maximum likelihood estimator, or MLE, is any maximizer of this likelihood:

$$
\hat\theta_n\in\operatorname*{argmax}_{\theta\in\Theta}L_n(\theta).
$$

Since the logarithm is strictly increasing, maximizing $L_n(\theta)$ is equivalent to maximizing the log-likelihood

$$
\ell_n(\theta)=\log L_n(\theta).
$$

Hence

$$
\hat\theta_n\in\operatorname*{argmax}_{\theta\in\Theta}\ell_n(\theta).
$$

If $X_1,\dots,X_n$ are i.i.d. with density $f_\theta$, then

$$
L_n(\theta)=\prod_{i=1}^n f_\theta(X_i),
$$

and therefore

$$
\ell_n(\theta)=\sum_{i=1}^n\log f_\theta(X_i).
$$

If $\ell_n(\theta)$ is differentiable and the maximizer is an interior point of $\Theta$, then the MLE usually satisfies the score equation

$$
\frac{\partial \ell_n(\theta)}{\partial\theta}=0.
$$

For vector-valued $\theta$, this becomes

$$
\nabla_\theta \ell_n(\theta)=0.
$$

The score equation is a necessary first-order condition for an interior maximum. It is not by itself sufficient; boundary points and multiple local maxima must still be checked.

## Invariance under transformation

Suppose $\theta\in\Theta$ is the parameter of interest. The MLE of $\theta$ is $\hat\theta=\text{argmax}_{\theta\in\Theta}p(x\mid\theta)$.

Consider $\eta=g(\theta)$, where $g(\cdot)$ is a certain transformation, it's easy to show that MLE of $\eta$ is given by $g(\hat\theta)$.

## Asymptotic Setup

Assume that

$$
X_1,\dots,X_n\overset{iid}{\sim}P_{\theta_0},
$$

where $\theta_0$ is the true parameter. The log-likelihood is

$$
\ell_n(\theta)=\sum_{i=1}^n\log f_\theta(X_i),
$$

and it is often useful to work with the average log-likelihood

$$
M_n(\theta)=\frac{1}{n}\ell_n(\theta)
=\frac{1}{n}\sum_{i=1}^n\log f_\theta(X_i).
$$

Its population version is

$$
M(\theta)=\mathbb E_{\theta_0}\log f_\theta(X).
$$

The MLE can equivalently be written as

$$
\hat\theta_n\in\operatorname*{argmax}_{\theta\in\Theta}M_n(\theta).
$$

## Consistency

### Theorem

Suppose the following conditions hold:

1. $\Theta$ is compact.
2. $\theta_0$ is identifiable, meaning

$$
P_\theta=P_{\theta_0}\quad\Longrightarrow\quad \theta=\theta_0.
$$

3. For almost every $x$, $\theta\mapsto \log f_\theta(x)$ is continuous on $\Theta$.
4. There exists an integrable envelope $M(x)$ such that

$$
\sup_{\theta\in\Theta}\left|\log f_\theta(x)\right|\le M(x),
\qquad
\mathbb E_{\theta_0}M(X)<\infty.
$$

Then

$$
\hat\theta_n\overset{p}\longrightarrow \theta_0.
$$

<u>Proof</u>

The proof has three steps:

1. Show that $M(\theta)$ is continuous.
2. Show that $M(\theta)$ has a unique maximizer at $\theta_0$.
3. Use uniform convergence of $M_n(\theta)$ to $M(\theta)$ to transfer the maximizer from the population criterion to the sample criterion.

First, since $\theta\mapsto\log f_\theta(x)$ is continuous and

$$
\sup_{\theta\in\Theta}\left|\log f_\theta(X)\right|\le M(X),
\qquad
\mathbb E_{\theta_0}M(X)<\infty,
$$

the dominated convergence theorem implies that

$$
M(\theta)=\mathbb E_{\theta_0}\log f_\theta(X)
$$

is continuous in $\theta$.

Next compare $M(\theta)$ with $M(\theta_0)$. We have

$$
M(\theta)-M(\theta_0)
=
\mathbb E_{\theta_0}\log \frac{f_\theta(X)}{f_{\theta_0}(X)}.
$$

By Jensen's inequality,

$$
\mathbb E_{\theta_0}\log \frac{f_\theta(X)}{f_{\theta_0}(X)}
\le
\log\mathbb E_{\theta_0}\frac{f_\theta(X)}{f_{\theta_0}(X)}.
$$

Since the expectation is taken under $P_{\theta_0}$,

$$
\mathbb E_{\theta_0}\frac{f_\theta(X)}{f_{\theta_0}(X)}
=
\int_{\{f_{\theta_0}>0\}}\frac{f_\theta(x)}{f_{\theta_0}(x)}f_{\theta_0}(x)\,d\mu(x)
=
\int_{\{f_{\theta_0}>0\}} f_\theta(x)\,d\mu(x)
\le 1.
$$

Therefore

$$
M(\theta)-M(\theta_0)\le 0,
$$

so

$$
M(\theta)\le M(\theta_0).
$$

Thus $\theta_0$ is a maximizer of the population criterion.

Now consider equality. Equality in Jensen's inequality holds only when

$$
\frac{f_\theta(X)}{f_{\theta_0}(X)}
$$

is constant $P_{\theta_0}$-almost surely. Because its expectation is at most $1$, equality in the whole chain above forces this constant to be $1$ and also forces $f_\theta$ to put no mass outside the support of $f_{\theta_0}$. Hence

$$
f_\theta=f_{\theta_0}
\quad P_{\theta_0}\text{-almost surely},
$$

and therefore

$$
P_\theta=P_{\theta_0}.
$$

By identifiability, this implies

$$
\theta=\theta_0.
$$

So $\theta_0$ is the unique maximizer of $M(\theta)$.

Next, by the uniform law of large numbers stated in the appendix, the compactness, continuity, and domination assumptions imply

$$
\sup_{\theta\in\Theta}|M_n(\theta)-M(\theta)|
\overset{p}\longrightarrow 0.
$$

Fix $\varepsilon>0$ and define the closed set

$$
A_\varepsilon=\{\theta\in\Theta:d(\theta,\theta_0)\ge \varepsilon\}.
$$

Since $\Theta$ is compact and $A_\varepsilon$ is closed, $A_\varepsilon$ is compact. Since $M$ is continuous and has the unique maximizer $\theta_0$,

$$
\sup_{\theta\in A_\varepsilon}M(\theta)<M(\theta_0).
$$

Define the positive separation

$$
\eta_\varepsilon
=
M(\theta_0)-\sup_{\theta\in A_\varepsilon}M(\theta)
>0.
$$

On the event

$$
\sup_{\theta\in\Theta}|M_n(\theta)-M(\theta)|<\frac{\eta_\varepsilon}{2},
$$

for every $\theta\in A_\varepsilon$,

$$
M_n(\theta)
\le
M(\theta)+\frac{\eta_\varepsilon}{2}
\le
M(\theta_0)-\eta_\varepsilon+\frac{\eta_\varepsilon}{2}
=
M(\theta_0)-\frac{\eta_\varepsilon}{2},
$$

while

$$
M_n(\theta_0)
\ge
M(\theta_0)-\frac{\eta_\varepsilon}{2}.
$$

Thus no maximizer of $M_n$ can lie in $A_\varepsilon$. Therefore

$$
\mathbb P_{\theta_0}\left(d(\hat\theta_n,\theta_0)\ge \varepsilon\right)
\le
\mathbb P_{\theta_0}\left(
\sup_{\theta\in\Theta}|M_n(\theta)-M(\theta)|
\ge
\frac{\eta_\varepsilon}{2}
\right)
\to 0.
$$

Hence, for every $\varepsilon>0$,

$$
\mathbb P_{\theta_0}\left(d(\hat\theta_n,\theta_0)\ge \varepsilon\right)\to 0.
$$

This proves

$$
\hat\theta_n\overset{p}\longrightarrow \theta_0.
$$

## Asymptotic Normality: One-Dimensional Case

Assume $\theta\in\Theta\subseteq\mathbb R$ and

$$
X_1,\dots,X_n\overset{iid}{\sim} f_{\theta_0}.
$$

We want to prove

$$
\sqrt n(\hat\theta_n-\theta_0)
\overset{d}\longrightarrow
N\left(0,\frac{1}{I(\theta_0)}\right),
$$

where $I(\theta_0)$ is the Fisher information for one observation.

### Regularity Conditions

A standard set of sufficient regularity conditions is:

1. $\theta_0\in\operatorname{int}(\Theta)$.
2. The MLE is consistent:

$$
\hat\theta_n\overset{p}\longrightarrow \theta_0.
$$

3. For $x$ in the sample space, $\log f_\theta(x)$ is three times continuously differentiable in $\theta$ in a neighborhood of $\theta_0$.
4. Differentiation and integration can be interchanged, so that

$$
\mathbb E_{\theta_0}\left[
\frac{\partial}{\partial\theta}\log f_{\theta_0}(X)
\right]=0.
$$

5. The Fisher information exists and satisfies

$$
I(\theta_0)
=
\mathbb E_{\theta_0}\left[
\left(
\frac{\partial}{\partial\theta}\log f_{\theta_0}(X)
\right)^2
\right]
=
-\mathbb E_{\theta_0}\left[
\frac{\partial^2}{\partial\theta^2}\log f_{\theta_0}(X)
\right],
$$

with

$$
0<I(\theta_0)<\infty.
$$

6. There exists an integrable function $B(X)$ such that, for $\theta$ near $\theta_0$,

$$
\left|
\frac{\partial^3}{\partial\theta^3}\log f_\theta(X)
\right|
\le B(X),
\qquad
\mathbb E_{\theta_0}B(X)<\infty.
$$

The last condition is mainly used to control the second derivative uniformly near $\theta_0$.

### Information Identities

Let the one-observation score be

$$
s_\theta(X)=\frac{\partial}{\partial\theta}\log f_\theta(X).
$$

Since

$$
\int f_\theta(x)\,d\mu(x)=1,
$$

differentiating under the integral sign gives

$$
0
=
\frac{\partial}{\partial\theta}\int f_\theta(x)\,d\mu(x)
=
\int \frac{\partial}{\partial\theta}f_\theta(x)\,d\mu(x).
$$

Because

$$
\frac{\partial}{\partial\theta}f_\theta(x)
=
f_\theta(x)\frac{\partial}{\partial\theta}\log f_\theta(x),
$$

we obtain

$$
\mathbb E_\theta s_\theta(X)=0.
$$

Differentiating once more,

$$
0
=
\int \frac{\partial^2}{\partial\theta^2}f_\theta(x)\,d\mu(x).
$$

Using

$$
\frac{\partial^2}{\partial\theta^2}f_\theta(x)
=
f_\theta(x)
\left[
\frac{\partial^2}{\partial\theta^2}\log f_\theta(x)
+
\left(
\frac{\partial}{\partial\theta}\log f_\theta(x)
\right)^2
\right],
$$

we get

$$
\mathbb E_\theta\left[
\frac{\partial^2}{\partial\theta^2}\log f_\theta(X)
\right]
+
\mathbb E_\theta\left[s_\theta(X)^2\right]
=0.
$$

Therefore

$$
I(\theta)
=
\mathbb E_\theta s_\theta(X)^2
=
-\mathbb E_\theta\left[
\frac{\partial^2}{\partial\theta^2}\log f_\theta(X)
\right].
$$

<u>Proof</u>

Since $\hat\theta_n$ is an interior maximizer, the first-order condition gives

$$
\ell_n'(\hat\theta_n)=0.
$$

Apply Taylor's theorem to $\ell_n'(\hat\theta_n)$ around $\theta_0$:

$$
0
=
\ell_n'(\hat\theta_n)
=
\ell_n'(\theta_0)
+
\ell_n''(\tilde\theta_n)(\hat\theta_n-\theta_0),
$$

where $\tilde\theta_n$ lies between $\hat\theta_n$ and $\theta_0$. Rearranging,

$$
\sqrt n(\hat\theta_n-\theta_0)
=
-
\frac{\frac{1}{\sqrt n}\ell_n'(\theta_0)}
{\frac{1}{n}\ell_n''(\tilde\theta_n)}.
$$

Now study the numerator and denominator separately.

For the numerator,

$$
\ell_n'(\theta_0)
=
\sum_{i=1}^n
\frac{\partial}{\partial\theta}\log f_{\theta_0}(X_i)
=
\sum_{i=1}^n s_{\theta_0}(X_i).
$$

The random variables $s_{\theta_0}(X_i)$ are i.i.d. with

$$
\mathbb E_{\theta_0}s_{\theta_0}(X)=0
$$

and

$$
\operatorname{Var}_{\theta_0}(s_{\theta_0}(X))=I(\theta_0).
$$

By the central limit theorem,

$$
\frac{1}{\sqrt n}\ell_n'(\theta_0)
=
\frac{1}{\sqrt n}\sum_{i=1}^n s_{\theta_0}(X_i)
\overset{d}\longrightarrow
N(0,I(\theta_0)).
$$

For the denominator,

$$
\frac{1}{n}\ell_n''(\tilde\theta_n)
=
\frac{1}{n}\sum_{i=1}^n
\frac{\partial^2}{\partial\theta^2}\log f_{\tilde\theta_n}(X_i).
$$

Since $\hat\theta_n\overset{p}\to\theta_0$ and $\tilde\theta_n$ lies between $\hat\theta_n$ and $\theta_0$,

$$
\tilde\theta_n\overset{p}\longrightarrow\theta_0.
$$

By the uniform law of large numbers applied to the class

$$
\left\{
\frac{\partial^2}{\partial\theta^2}\log f_\theta:\theta\text{ near }\theta_0
\right\},
$$

we have

$$
\sup_{\theta\text{ near }\theta_0}
\left|
\frac{1}{n}\sum_{i=1}^n
\frac{\partial^2}{\partial\theta^2}\log f_\theta(X_i)
-
\mathbb E_{\theta_0}
\frac{\partial^2}{\partial\theta^2}\log f_\theta(X)
\right|
\overset{p}\longrightarrow 0.
$$

Therefore,

$$
\frac{1}{n}\ell_n''(\tilde\theta_n)
-
\mathbb E_{\theta_0}
\left[
\frac{\partial^2}{\partial\theta^2}\log f_{\tilde\theta_n}(X)
\right]
\overset{p}\longrightarrow 0.
$$

By continuity and dominated convergence,

$$
\mathbb E_{\theta_0}
\left[
\frac{\partial^2}{\partial\theta^2}\log f_{\tilde\theta_n}(X)
\right]
\overset{p}\longrightarrow
\mathbb E_{\theta_0}
\left[
\frac{\partial^2}{\partial\theta^2}\log f_{\theta_0}(X)
\right].
$$

Using the information identity,

$$
\mathbb E_{\theta_0}
\left[
\frac{\partial^2}{\partial\theta^2}\log f_{\theta_0}(X)
\right]
=
-I(\theta_0).
$$

Thus

$$
\frac{1}{n}\ell_n''(\tilde\theta_n)
\overset{p}\longrightarrow
-I(\theta_0).
$$

Combining the numerator and denominator,

$$
\frac{1}{\sqrt n}\ell_n'(\theta_0)
\overset{d}\longrightarrow
N(0,I(\theta_0)),
$$

and

$$
\frac{1}{n}\ell_n''(\tilde\theta_n)
\overset{p}\longrightarrow
-I(\theta_0).
$$

By Slutsky's theorem,

$$
\sqrt n(\hat\theta_n-\theta_0)
=
-
\frac{\frac{1}{\sqrt n}\ell_n'(\theta_0)}
{\frac{1}{n}\ell_n''(\tilde\theta_n)}
\overset{d}\longrightarrow
\frac{Z}{I(\theta_0)},
$$

where

$$
Z\sim N(0,I(\theta_0)).
$$

Therefore

$$
\sqrt n(\hat\theta_n-\theta_0)
\overset{d}\longrightarrow
N\left(0,\frac{1}{I(\theta_0)}\right).
$$

## Asymptotic Normality: Multidimensional Case

Now let

$$
\theta\in\Theta\subseteq\mathbb R^d.
$$

The goal is

$$
\sqrt n(\hat\theta_n-\theta_0)
\overset{d}\longrightarrow
N_d\left(0,I(\theta_0)^{-1}\right),
$$

where $I(\theta_0)$ is the $d\times d$ Fisher information matrix.

### Notation

Define the score vector

$$
S_n(\theta)
=
\nabla_\theta \ell_n(\theta)
=
\sum_{i=1}^n\nabla_\theta\log f_\theta(X_i),
$$

and the Hessian matrix

$$
H_n(\theta)
=
\nabla_\theta^2 \ell_n(\theta)
=
\sum_{i=1}^n\nabla_\theta^2\log f_\theta(X_i).
$$

For one observation, write

$$
s_\theta(X)=\nabla_\theta\log f_\theta(X).
$$

The Fisher information matrix is

$$
I(\theta_0)
=
\mathbb E_{\theta_0}\left[
s_{\theta_0}(X)s_{\theta_0}(X)^\top
\right].
$$

Under the usual regularity conditions,

$$
\mathbb E_{\theta_0}s_{\theta_0}(X)=0,
$$

and

$$
I(\theta_0)
=
-\mathbb E_{\theta_0}\left[
\nabla_\theta^2\log f_{\theta_0}(X)
\right].
$$

Assume $I(\theta_0)$ is nonsingular.

<u>Proof</u>

Since $\hat\theta_n$ is an interior maximizer,

$$
S_n(\hat\theta_n)=0.
$$

By the multivariate mean value form of Taylor's theorem,

$$
0
=
S_n(\hat\theta_n)
=
S_n(\theta_0)
+
\bar H_n(\hat\theta_n-\theta_0),
$$

where

$$
\bar H_n
=
\int_0^1
H_n\left(\theta_0+t(\hat\theta_n-\theta_0)\right)
\,dt.
$$

Thus

$$
\sqrt n(\hat\theta_n-\theta_0)
=
-
\left(\frac{1}{n}\bar H_n\right)^{-1}
\left(\frac{1}{\sqrt n}S_n(\theta_0)\right),
$$

whenever $\frac{1}{n}\bar H_n$ is invertible. Since its probability limit is nonsingular, this invertibility holds with probability tending to $1$.

First consider the score term:

$$
\frac{1}{\sqrt n}S_n(\theta_0)
=
\frac{1}{\sqrt n}
\sum_{i=1}^n s_{\theta_0}(X_i).
$$

The vectors $s_{\theta_0}(X_i)$ are i.i.d. with mean zero and covariance matrix $I(\theta_0)$. By the multivariate central limit theorem,

$$
\frac{1}{\sqrt n}S_n(\theta_0)
\overset{d}\longrightarrow
N_d(0,I(\theta_0)).
$$

Next consider the Hessian term. Since

$$
\hat\theta_n\overset{p}\longrightarrow\theta_0,
$$

all points on the line segment between $\theta_0$ and $\hat\theta_n$ also converge in probability to $\theta_0$. By the uniform law of large numbers applied entrywise to the Hessian class,

$$
\sup_{\theta\in U}
\left\|
\frac{1}{n}H_n(\theta)
-
\mathbb E_{\theta_0}\nabla_\theta^2\log f_\theta(X)
\right\|
\overset{p}\longrightarrow 0
$$

for some neighborhood $U$ of $\theta_0$. Hence

$$
\frac{1}{n}\bar H_n
-
\int_0^1
\mathbb E_{\theta_0}\nabla_\theta^2
\log f_{\theta_0+t(\hat\theta_n-\theta_0)}(X)
\,dt
\overset{p}\longrightarrow 0.
$$

By continuity and dominated convergence,

$$
\int_0^1
\mathbb E_{\theta_0}\nabla_\theta^2
\log f_{\theta_0+t(\hat\theta_n-\theta_0)}(X)
\,dt
\overset{p}\longrightarrow
\mathbb E_{\theta_0}\nabla_\theta^2\log f_{\theta_0}(X).
$$

Using the information identity,

$$
\mathbb E_{\theta_0}\nabla_\theta^2\log f_{\theta_0}(X)
=
-I(\theta_0).
$$

Therefore

$$
\frac{1}{n}\bar H_n
\overset{p}\longrightarrow
-I(\theta_0).
$$

Since $I(\theta_0)$ is nonsingular, the continuous mapping theorem gives

$$
\left(\frac{1}{n}\bar H_n\right)^{-1}
\overset{p}\longrightarrow
\left[-I(\theta_0)\right]^{-1}
=
-I(\theta_0)^{-1}.
$$

Consequently,

$$
-
\left(\frac{1}{n}\bar H_n\right)^{-1}
\overset{p}\longrightarrow
I(\theta_0)^{-1}.
$$

Combining this with the multivariate CLT and applying Slutsky's theorem,

$$
\sqrt n(\hat\theta_n-\theta_0)
\overset{d}\longrightarrow
I(\theta_0)^{-1}Z,
$$

where

$$
Z\sim N_d(0,I(\theta_0)).
$$

Since a linear transformation of a multivariate normal vector is still multivariate normal,

$$
I(\theta_0)^{-1}Z
\sim
N_d\left(
0,
I(\theta_0)^{-1}I(\theta_0)I(\theta_0)^{-1}
\right).
$$

Thus

$$
\sqrt n(\hat\theta_n-\theta_0)
\overset{d}\longrightarrow
N_d\left(0,I(\theta_0)^{-1}\right).
$$

## Likelihood ratio test

The likelihood method also provides a natural way to construct hypothesis tests.

Suppose we want to test

$$
H_0:\theta\in\Theta_0
\quad\text{vs.}\quad
H_1:\theta\in\Theta.
$$

We construct the likelihood statistic as follows.

Let

$$
\Lambda_n
=
\frac{\sup_{\theta\in\Theta_0}L_n(\theta)}
{\sup_{\theta\in\Theta}L_n(\theta)}.
$$

If the true parameter $\theta$ lies in $\Theta_0\subset\Theta$, then the numerator and denominator should be the same, $\Lambda_n$ close to $1$ accordingly. If the true parameter $\theta$ lies out of $\Theta_0$, then the numerator should be much less the denominator, making $\Lambda_n$ small.

Actually, under certain regularity conditions, we have Wilks' theorem.

## Wilks' theorem for the likelihood ratio test

Let $X_1,\dots,X_n$ be i.i.d. from a regular parametric family $\{P_\theta:\theta\in\Theta\subset\mathbb R^p\}$. Suppose we test

$$
H_0:\theta\in\Theta_0
$$

where $\Theta_0$ is a smooth submanifold of $\Theta$ with dimension $q<p$. Let

$$
\Lambda_n
=
\frac{\sup_{\theta\in\Theta_0} L_n(\theta)}
{\sup_{\theta\in\Theta} L_n(\theta)}.
$$

If the true parameter $\theta_0\in\Theta_0$, and standard regularity conditions hold, then

$$
-2\log\Lambda_n
\overset{d}\longrightarrow
\chi^2_{p-q}.
$$

The degrees of freedom equal the number of restrictions imposed by $H_0$.





## Appendix

### ULLN

The consistency proof uses the following uniform law of large numbers.

Let $\Theta$ be compact. Suppose $g(x,\theta)$ satisfies:

1. For almost every $x$, $\theta\mapsto g(x,\theta)$ is continuous on $\Theta$.
2. There exists an integrable envelope $G(x)$ such that

$$
\sup_{\theta\in\Theta}|g(x,\theta)|\le G(x),
\qquad
\mathbb EG(X)<\infty.
$$

Define

$$
P_ng_\theta=\frac{1}{n}\sum_{i=1}^n g(X_i,\theta),
\qquad
Pg_\theta=\mathbb E g(X,\theta).
$$

Then

$$
\sup_{\theta\in\Theta}|P_ng_\theta-Pg_\theta|
\overset{p}\longrightarrow 0.
$$

The same conclusion also holds almost surely under the same conditions.

<u>Proof</u>

The main idea is that compactness lets us approximate the whole parameter space by finitely many representative points, and domination makes the approximation valid in expectation.

Fix $\varepsilon>0$. For a fixed $\theta_0\in\Theta$, define

$$
h_\delta(X)
=
\sup_{\theta\in\Theta:\|\theta-\theta_0\|\le\delta}
|g(X,\theta)-g(X,\theta_0)|.
$$

Because $g(X,\theta)$ is continuous in $\theta$ for almost every $X$,

$$
h_\delta(X)\downarrow 0
\quad\text{as}\quad
\delta\downarrow 0
$$

for almost every $X$. Also,

$$
0\le h_\delta(X)\le 2G(X),
$$

and $2G(X)$ is integrable. Hence, by dominated convergence,

$$
\mathbb E h_\delta(X)\to 0.
$$

Therefore, for each $\theta_0\in\Theta$, there exists a neighborhood $U_{\theta_0}$ such that

$$
\mathbb E\left[
\sup_{\theta\in U_{\theta_0}}
|g(X,\theta)-g(X,\theta_0)|
\right]
<
\varepsilon.
$$

The collection $\{U_{\theta_0}:\theta_0\in\Theta\}$ is an open cover of $\Theta$. Since $\Theta$ is compact, there exist finitely many points $\theta_1,\dots,\theta_m$ such that

$$
\Theta\subseteq U_{\theta_1}\cup\cdots\cup U_{\theta_m}.
$$

For each $j$, define the local oscillation

$$
H_j(X)
=
\sup_{\theta\in U_{\theta_j}}
|g(X,\theta)-g(X,\theta_j)|.
$$

By construction,

$$
\mathbb EH_j(X)<\varepsilon.
$$

For every $\theta\in\Theta$, choose $j(\theta)$ such that $\theta\in U_{\theta_{j(\theta)}}$. Then

$$
|P_ng_\theta-Pg_\theta|
\le
|P_ng_{\theta_{j(\theta)}}-Pg_{\theta_{j(\theta)}}|
+
P_nH_{j(\theta)}
+
PH_{j(\theta)}.
$$

Taking the supremum over $\theta$ gives

$$
\sup_{\theta\in\Theta}|P_ng_\theta-Pg_\theta|
\le
\max_{1\le j\le m}|P_ng_{\theta_j}-Pg_{\theta_j}|
+
\max_{1\le j\le m}P_nH_j
+
\max_{1\le j\le m}PH_j.
$$

The first term converges to $0$ by the ordinary law of large numbers applied to the finite collection $g_{\theta_1},\dots,g_{\theta_m}$. For the second term, the ordinary law of large numbers gives

$$
P_nH_j\to PH_j
$$

for each fixed $j$, and therefore

$$
\max_{1\le j\le m}P_nH_j
\to
\max_{1\le j\le m}PH_j
<
\varepsilon.
$$

The third term is also less than $\varepsilon$. Thus

$$
\limsup_{n\to\infty}
\sup_{\theta\in\Theta}|P_ng_\theta-Pg_\theta|
\le 2\varepsilon
$$

in probability, and in fact almost surely. Since $\varepsilon>0$ is arbitrary,

$$
\sup_{\theta\in\Theta}|P_ng_\theta-Pg_\theta|
\to 0.
$$

This proves the uniform law of large numbers.

### Proof of Wilks' theorem

Write the log-likelihood as

$$
\ell_n(\theta)=\sum_{i=1}^n \log f_\theta(X_i).
$$

Let $\hat\theta_n$ be the unrestricted MLE and $\tilde\theta_n$ the restricted MLE under $H_0$. Then

$$
-2\log\Lambda_n
=
2\{\ell_n(\hat\theta_n)-\ell_n(\tilde\theta_n)\}.
$$

The proof is local. Since $\Theta_0$ is a smooth $q$-dimensional submanifold, we may use a smooth reparametrization near $\theta_0$ so that

$$
\theta=(\alpha,\beta),
\qquad
\alpha\in\mathbb R^q,
\qquad
\beta\in\mathbb R^{p-q},
$$

and the null hypothesis becomes

$$
H_0:\theta_{q+1}=\cdots=\theta_p=0.
$$

Equivalently, $H_0$ is $\beta=0$. The true value has the form

$$
\theta_0=(\alpha_0,0).
$$

Let

$$
S_n(\theta_0)=\nabla_\theta \ell_n(\theta_0),
\qquad
J_n(\theta_0)=-\nabla_\theta^2\ell_n(\theta_0).
$$

Here $S_n(\theta_0)$ is the score and $J_n(\theta_0)$ is the observed information. Under the usual regularity conditions,

$$
\Delta_n
:=
\frac{1}{\sqrt n}S_n(\theta_0)
\overset{d}\longrightarrow
\Delta\sim N_p(0,I),
$$

and

$$
\frac{1}{n}J_n(\theta_0)
\overset{p}\longrightarrow
I,
$$

where $I=I(\theta_0)$ is the Fisher information matrix.

Now look at local alternatives of the form

$$
\theta=\theta_0+\frac{h}{\sqrt n},
\qquad
h\in\mathbb R^p.
$$

A second-order Taylor expansion gives

$$
\ell_n\left(\theta_0+\frac{h}{\sqrt n}\right)-\ell_n(\theta_0)
=
h^\top\Delta_n
-\frac{1}{2}h^\top
\left\{\frac{1}{n}J_n(\theta_0)\right\}
h
+r_n(h).
$$

The remainder satisfies

$$
\sup_{\|h\|\le C}|r_n(h)|
\overset{p}\longrightarrow 0
$$

for every fixed $C<\infty$. This follows from the uniform convergence of the Hessian in a shrinking neighborhood of $\theta_0$.

Therefore, on every bounded set of $h$, the local log-likelihood is asymptotically

$$
Q_n(h)
=
h^\top\Delta_n-\frac{1}{2}h^\top I h+o_p(1),
$$

and the limiting quadratic criterion is

$$
Q(h)
=
h^\top\Delta-\frac{1}{2}h^\top I h.
$$

The unrestricted local parameter is

$$
\hat h_n=\sqrt n(\hat\theta_n-\theta_0),
$$

and the restricted local parameter is

$$
\tilde h_n=\sqrt n(\tilde\theta_n-\theta_0).
$$

Both are $O_p(1)$ under the regularity conditions and under $H_0$, so the local quadratic expansion applies to them. Thus the likelihood ratio statistic has the same limit as

$$
2\left\{
\sup_{h\in\mathbb R^p}Q(h)
-
\sup_{h\in T_0}Q(h)
\right\},
$$

where the tangent space of the null hypothesis is

$$
T_0
=
\{h\in\mathbb R^p:h_{q+1}=\cdots=h_p=0\}.
$$

The unrestricted maximizer of $Q(h)$ solves

$$
\nabla_h Q(h)=\Delta-Ih=0,
$$

so

$$
\hat h=I^{-1}\Delta.
$$

The restricted maximizer is the projection of $\hat h$ onto $T_0$ using the Fisher-information inner product

$$
\langle a,b\rangle_I=a^\top I b.
$$

This means that

$$
2\left\{
Q(\hat h)-Q(\tilde h)
\right\}
=
(\hat h-\tilde h)^\top I(\hat h-\tilde h),
$$

the squared distance from the unrestricted limit to the null tangent space in the Fisher-information metric.

To see the distribution of this squared distance, complete the square:

$$
Q(h)
=
\frac{1}{2}\Delta^\top I^{-1}\Delta
-
\frac{1}{2}
(h-I^{-1}\Delta)^\top I(h-I^{-1}\Delta).
$$

Thus maximizing $Q(h)$ is the same as minimizing the Fisher-information distance from $h$ to $I^{-1}\Delta$.

Now whiten the Fisher information. Let

$$
y=I^{1/2}h,
\qquad
W=I^{-1/2}\Delta.
$$

Since $\Delta\sim N_p(0,I)$,

$$
W\sim N_p(0,I_p).
$$

The null tangent space $T_0$ is transformed into the $q$-dimensional linear subspace

$$
V=I^{1/2}T_0.
$$

In these coordinates,

$$
2\{Q(\hat h)-Q(\tilde h)\}
=
\operatorname{dist}(W,V)^2,
$$

where $\operatorname{dist}(W,V)$ is the ordinary Euclidean distance from $W$ to the subspace $V$.

Finally, because $W$ is standard normal and therefore rotationally invariant, we may apply an orthogonal rotation that maps $V$ to the coordinate subspace spanned by the first $q$ basis vectors. This rotation does not change the distribution of $W$. In the rotated coordinates, write $W=(W_1,\dots,W_p)^\top$. Then

$$
\operatorname{dist}(W,V)^2
=
W_{q+1}^2+\cdots+W_p^2.
$$

The random variables $W_{q+1},\dots,W_p$ are independent standard normal variables. Hence

$$
W_{q+1}^2+\cdots+W_p^2
\sim
\chi^2_{p-q}.
$$

Consequently,

$$
-2\log\Lambda_n
\overset{d}\longrightarrow
\chi^2_{p-q}.
$$

This proves Wilks' theorem.
