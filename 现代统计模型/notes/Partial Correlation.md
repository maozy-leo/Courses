Partial Correlation



Correlation is defined to describe the linear relationship between two variables. And when we have more than two R.V., the relationship between variables become more complex. Thus we introduce the following definition to describe the conditioned linear relationship between R.V.

Formally, 
$$
X=(X_1,\dots,X_p)',\,\,\,\, p\ge3
$$
is a random vector with finite second-order moment $Var(X)=\Sigma$. Write $Z=(X_3,\dots,X_p)', A=(X_1,X_2)'$.

<u>Def1</u>

Let $\Omega=\Sigma^{-1}=(\omega_{ij})$ be the precision matrix. Then the partial correlation between $X_1,X_2$ conditioned on $X_3,\dots,X_p$ is defined as
$$
\rho_{12,-12}=-\dfrac{\omega_{12}}{\sqrt{\omega_{11}\omega_{22}}}
$$
<u>Def2</u>

Partition $\Sigma$ by $A,Z$ into
$$
\Sigma=\begin{pmatrix}
\Sigma_{AA}&\Sigma_{AZ}\\
\Sigma_{ZA}&\Sigma_{ZZ}
\end{pmatrix}
$$
Define the Schur-complement matrix as
$$
S=\Sigma_{AA}-\Sigma_{AZ}\Sigma_{ZZ}^{-1}\Sigma_{ZA}
$$
Write $S$ as
$$
S=\begin{pmatrix}
s_{11}&s_{12}\\
s_{21}&s_{22}
\end{pmatrix}
$$
Then the partial correlation is defined as
$$
\rho_{12,-12}=\dfrac{s_{12}}{\sqrt{s_{11}s_{22}}}
$$
<u>Def3</u>

Project $X_1,X_2$ linearly into the span of $Z$, denoted as $\hat X_{1},\hat X_2$.

Define the residuals as
$$
R_1=X_1-\hat X_1, R_2=X_2-\hat X_2
$$
Then the partial correlation of $X_1,X_2$ is defined as 
$$
\rho_{12,-12}=corr(R_1,R_2)=\dfrac{Cov(R_1,R_2)}{\sqrt{Var(R_1)Var(R_2)}}
$$


We now prove the equivalence of the above definitions.

<u>Proof</u>

We first find the explicit form of $\hat X_1,\hat X_2$.
$$
\begin{align}
\hat X_1&=\text{argmin}_{\tilde Z\in\text{span}Z}\mathbb E(X_1-\tilde Z)^2\\
&=[\text{argmin}_{\alpha\in\R^{p-2}}\mathbb E(X_1-\alpha' Z)^2]'Z\\
&=\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}Z

\end{align}
$$
Similarly,
$$
\hat X_2=\Sigma_{X_2Z}\Sigma_{ZZ}^{-1}Z
$$
Thus,
$$
\begin{align}
Cov(R_1,R_2)&=Cov(X_1-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}Z,X_2-\Sigma_{X_2Z}\Sigma_{ZZ}^{-1}Z)\\
&=\sigma_{12}-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}\Sigma_{ZX_2}

\end{align}
$$

$$
\begin{align}
Var(R_1)&=Var(X_1-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}Z)\\
&=\sigma_1^2-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}\Sigma_{ZX_1}

\end{align}
$$

Similarly,
$$
Var(R_2)=\sigma^2_2-\Sigma_{X_2Z}\Sigma^{-1}_{ZZ}\Sigma_{ZX_2}
$$
From the inverse formula of blocked matrix, the upper-left block of $\Omega$ is the inverse of $\Sigma_{AA}-\Sigma_{AZ}\Sigma^{-1}_{ZZ}\Sigma_{ZA}$.

In detail,  $\Sigma_{AA}-\Sigma_{AZ}\Sigma^{-1}_{ZZ}\Sigma_{ZA}$ can be written as
$$
\begin{pmatrix}
w_{11}&w_{12}\\
w_{21}&w_{22}
\end{pmatrix}
=
\begin{pmatrix}
\sigma^2_1-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}\Sigma_{ZX_1}&
\sigma_{12}-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}\Sigma_{ZX_2}\\
\sigma_{12}-\Sigma_{X_1Z}\Sigma_{ZZ}^{-1}\Sigma_{ZX_2}&
\sigma^2_2-\Sigma_{X_2Z}\Sigma_{ZZ}^{-1}\Sigma_{ZX_2}
\end{pmatrix}
=
\begin{pmatrix}
Var(R_1)&Cov(R_1,R_2)\\
Cov(R_1,R_2)&Var(R_2)
\end{pmatrix}
$$
Take inverse, we have
$$
\begin{pmatrix}
w_{11}&w_{12}\\
w_{21}&w_{22}
\end{pmatrix}^{-1}=
\dfrac{1}{\det S}\begin{pmatrix}
w_{22}&-w_{21}\\
-w_{12}&w_{11}
\end{pmatrix}
$$
Thus,
$$
\begin{align}
\rho_{12,-12}&=-\dfrac{w_{12}}{\sqrt{w_{11}w_{22}}}\\
&=\dfrac{Cov(R_1,R_2)/\det S}{\sqrt{ Var(R_1)/\det S\cdot Var(R_2)/\det S  }}  \\
&=\dfrac{Cov(R_1,R_2)}{\sqrt{Var(R_1)Var(R_2)}}
\end{align}
$$


In some specific cases, the partial correlation also have other forms of definition.

- $p=3$

In this case, $\Sigma_{ZZ}$ degenerates into a number. Thus
$$
Var(R_1)=\sigma_1^2-\dfrac{\sigma^2_{13}}{\sigma^2_{3}}\\
Var(R_2)=\sigma_2^2-\dfrac{\sigma_{23}^2}{\sigma_{3}^2}\\
Cov(R_1,R_2)=\sigma_{12}-\dfrac{\sigma_{13}\sigma_{23}}{\sigma_{3}^2}
$$
Thus, 
$$
\begin{align}
\rho_{12,-12}&=\dfrac{Cov(R_1,R_2)}{\sqrt{Var(R_1)Var(R_2)}}\\
&=\dfrac{\sigma_{12}\sigma_{3}^2-\sigma_{13}\sigma_{23}}{\sqrt{(\sigma_{13}^2-\sigma_{3}^2\sigma_{1}^2)(\sigma_{23}^2-\sigma_{3}^2\sigma_2^2)}}\\
&=\dfrac{\dfrac{\sigma_{12}}{\sigma_1\sigma_2}- \dfrac{\sigma_{13}}{\sigma_1\sigma_3} \dfrac{\sigma_{23}}{\sigma_2\sigma_3} }{\sqrt{( 1 -\dfrac{\sigma_{13}^2}{\sigma_1^2\sigma_{3}^2})(1-\dfrac{\sigma_{23}^2}{\sigma_2^2\sigma_{3}^2})}}\\
&=\dfrac{\rho_{12}-\rho_{13}\rho_{23}}{\sqrt{(1-\rho_{13}^2)(1-\rho_{23}^2)}}
\end{align}
$$
