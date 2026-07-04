# BLUP--best linear unbiased predictor

Consider linear model with random effects
$$
y=X\beta+Zu+\epsilon
$$
where
$$
\begin{pmatrix}
u\\\epsilon
\end{pmatrix}\sim
\mathcal N\left(
\begin{pmatrix}
0\\0
\end{pmatrix},
\begin{pmatrix}
G&0\\0&R
\end{pmatrix}

\right)
$$
Our goal is to find a predictor of the random effects $u$, given data $y$. We denote our predictor as $\hat u$.

Best linear unbiased predictor, in short BLUP, gives a predictor for the aforementioned question, which satisfies the following requirements:

- $\hat u$ is a linear function of $y$;
- $\hat u$ is unbiased for $u$, here the expectation is unconditional;
- $\text{Var}(\hat u-u)$ is no larger than $\text{Var}(\hat w-u)$, where $\hat w$ is any unbiased linear function of $u$.

We now derive the concrete form of BLUP $\hat u$ given the requirements.

First, by the decomposition of variance,
$$
\text{Var}(\hat u-u)=\mathbb E\text{Var}(\hat u-u\mid y)+\text{Var}\mathbb E(\hat u-u\mid y)
$$
 Thus, taken $\hat u=\mathbb E(u\mid y)$ yields the least variance, whose exact form will be derived next.

Notice $u, \epsilon$ are jointly Gaussian distributed. Thus $(u,y)^\top$ is Gaussian because it can be written as linear function of $u,\epsilon$.

We have
$$
\begin{pmatrix}
u\\y
\end{pmatrix}\sim
\mathcal N\left(
\begin{pmatrix}
0\\X\beta
\end{pmatrix},
\begin{pmatrix}
G&GZ^\top\\ZG&V
\end{pmatrix}

\right)
$$
where $V=ZGZ^\top+R$.

Then the Conditioned distribution of $u$ given $y$ is
$$
u\mid y\sim \mathcal N(GZ^{\top}V^{-1}(y-X\beta),G-ZGV^{-1}GZ^\top)
$$
Thus, BLUP is given by
$$
\hat u = GZ^\top V^{-1}(y-X\beta)
$$
In practice, we replace $\beta$ by $\hat\beta_{V}$ to approximate the real $\hat u$, which is called EBLUP. Here $\hat\beta_V$ is the GLS estimator
$$
\hat\beta_V=(X'V^{-1}X)^{-1}XV^{-1}y
$$


