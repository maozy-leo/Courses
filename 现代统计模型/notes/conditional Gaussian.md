In this note, we will derive the conditional distribution of multi-normal variables.

Suppose $X=(X_1,\dots,X_p)'\sim N(\mu,\Sigma)$ is a $p$ dimensional Gaussion vector.

We partition $X$ and its expectation/covariance matrix into the following form:
$$
X=\begin{pmatrix}
X_{(1)}\\X_{(2)}
\end{pmatrix}
\,\,\,\,\mu=\begin{pmatrix}
\mu_1\\\mu_2
\end{pmatrix}
\,\,\,\,\Sigma=\begin{pmatrix}
\Sigma_{11}&\Sigma_{12}\\
\Sigma_{21}&\Sigma_{22}
\end{pmatrix}
$$
with the first block row $q$, second row $p-q$ columns.

In the following part, we will caculate marginal conditional distribution $X_{(1)}|X_{(2)}=x_2$ by projection.



Suppose $X_{(1)}$ is in the following form, where $\varepsilon$ is the residual vector.
$$
X_{(1)}-\mu_1=B(X_{(2)}-\mu_2)+\varepsilon
$$
We have, 
$$
Cov(X_{(2)},\varepsilon)=0
$$
Then, expand $\varepsilon$
$$
0=Cov(X_{(2)},X_{(1)}-BX_{(2)})=Cov(X_{(2)},X_{(1)})-Cov(X_{(2)},X_{(2)})B'
$$
We can then get the solution of $B=\Sigma_{12}\Sigma_{22}^{-1}$.

$\begin{pmatrix}\varepsilon\\X_{(2)}\end{pmatrix}$ is a linear transformation of $X$, thus Guassian.

Their 0 covariance implies independence.

Thus, conditional distribution of $\varepsilon$ given $X_{(2)}$ is equal to its unconditional distribution.

From the definition of $\varepsilon$, we know that given $X_{(2)}$, $X_{(1)}$ can be generated from linear transformation of $\varepsilon$, thus Gaussian.

We only need to find its conditional expectation and variance.
$$
\begin{align}
\mathbb E(X_{(1)}|X_{(2)}=x_2)&=\mathbb E(\mu_1+B(X_{(2)}-\mu_2)+\varepsilon|X_{(2)}=x_2)\\
&=\mu_1+\Sigma_{12}\Sigma_{22}^{-1}(x_2-\mu_2)+\mathbb E(\varepsilon)\\
&=\mu_1+\Sigma_{12}\Sigma_{22}^{-1}(x_2-\mu_2)

\end{align}
$$

$$
\begin{align}
Var(X_{(1)}|X_2=x_2)&=Var(\mu_1+B(X_{(2)}-\mu_2)+\varepsilon|X_{(2)}=x_2)\\
&=Var(\varepsilon|X_{2}=x_2)\\
&=Var(\varepsilon)\\
&=Var(X_{(1)}-BX_{(2)})\\
&=\Sigma_{11}-2\Sigma_{12}B'+B\Sigma_{22}B'\\
&=\Sigma_{11}-\Sigma_{12}\Sigma_{22}^{-1}\Sigma_{12}

\end{align}
$$

Finally, we conclud that conditional on $X_{(2)}=x_2$, $X_{(1)}\sim N(\mu_1+\Sigma_{12}\Sigma_{22}^{-1}(x_2-\mu_2),\Sigma_{11}-\Sigma_{12}\Sigma_{22}^{-1}\Sigma_{12})$.