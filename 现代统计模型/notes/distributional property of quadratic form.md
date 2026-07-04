In this part of note, we mainly focus on the proof of the following two theorems:

Given a random vector $Y\sim N(\mu,\Sigma)$, where $\Sigma$ is a positive definite matrix.

<u>Thm 1</u>

$A$ is a $n\times n$ symmetric matrix. If $A\Sigma$ is idempotent, then we have $Y'AY\sim \chi^2_{m}(\mu'A\mu)$, where $\rank(A)=m$.

<u>Thm 2</u>

$A_1, A_2$ are two $n\times n$ symmetric matrixes. If $A_1\Sigma A_2=0$, then $Y'A_1Y\perp Y'A_2Y$.

<u>Thm 3</u>

$A$ is an $m\times n$ matrix, while $B$ is an $n\times n$ symmetric matrix. If $A\Sigma B=0$, then $AY\perp Y'BY$.

Notice: Thm1 and Thm2 can be strengthened to iff condition.



To complete the proofs, we first introduce two lemmas and definition of non-centered chi-square distribution.

<u>Lemma 1</u> *(Chi-square)*

Given Gaussian vector $X\sim N(\mu, I_n)$. $A$ is an $n\times n$ symmetric matrix. Then $X'AX\sim\chi_r(\mu'A\mu)$ iff $A$ is idempotent with $\rank$ $r$.

<u>Lemma 2</u> *(Craig-Sakamoto)*

Given Gaussian vector $X\sim N(\mu, I_n)$. $A,B$ are two $n\times n$ symmetrix matrixes. Then $X'AX\perp X'BX$ iff $AB=0$.

<u>Def</u> *(Non-contered chi-square distribution)*

Given independent Gaussian vectors $X_i\sim N(\mu_i,1),i=1,\dots,n$. Then we the squared sum $\sum_{i=1}^nX_i^2$ follows the non-centered chi-square distribution $\chi^2_n(\sum_{i=1}^n\mu_i^2)$ with degree of freedom $n$.

A random variable follows non-centered chi-square distribution $\chi^2_n(\lambda)$ iff its moment generating function has the following form:
$$
M(t)=\mathbb Ee^{tX}=(1-2t)^{-n/2}exp\{ \dfrac{\lambda t}{1-2t} \}
$$



We first prove the thm1, thm2 using the two lemmas. And the lemmas will be proved in later part.

<u>Proof of Thm 1</u>

Consider $X=\Sigma^{-1/2}Y\sim N(\Sigma^{-1/2}\mu, I_n)$.

Then $Y'AY=X'\Sigma^{1/2}A\Sigma^{1/2} X$. From Lemma1, we only need to check

- $\rank(\Sigma^{1/2}A\Sigma^{1/2})=m$, which is trivial
- $\Sigma^{1/2}A\Sigma^{1/2}$ is idempotent:

$(\Sigma^{1/2}A\Sigma^{1/2})^2=\Sigma^{1/2}A\Sigma^{1/2}$ $\Leftarrow$$\Sigma^{1/2}A\Sigma A\Sigma^{1/2}=\Sigma^{1/2}A\Sigma^{1/2}$$\Leftarrow$$A\Sigma$ is idempotent.

- the parameter of non-centered chi-square distribution is $\mu'A\mu$, which can be directly proved



<u>Proof of Thm 2</u>

Similarly consider $X=\Sigma^{-1/2}Y\sim N(\Sigma^{-1/2}\mu, I_n)$.

Then $Y'A_iY=X'\Sigma^{1/2}A_i\Sigma^{1/2} X, i=1,2$. From Lemma2, we only need to check

$\Sigma^{1/2}A_1\Sigma A_2\Sigma^{1/2}=0$, which directly comed from $A_1\Sigma A_2=0$



<u>Proof of Lemma 1</u>

Due to $A$ is a real symmetric matrix, we can decompose it by eigenvalues in the following form:
$$
A=U'\Lambda U=\sum_{i=1}^n\lambda_i u_iu_i'
$$
Thus,
$$
X'AX=\sum_{i=1}^n\lambda_iX'u_iu_i'X=\sum_{i=1}^n\lambda_i(u_i'X)'(u_i'X)
$$
We then caculate the moment generating function of $X'AX$.

Because $Cov(u_i'X,u_j'X)=0, i\ne j$, each term of $X'AX$ is independent.

Thus we can break the calculation of $X'AX$ down into the caculation of $\lambda_i(u_i'X)'(u_i'X)$:


$$
(1-2\lambda_it)^{-1/2}exp\{ \dfrac{\lambda_it}{1-2\lambda _it}(u_i'\mu)^2 \}
$$
Then, $X'AX$ has the following moment generating function:
$$
\begin{align}
M(t)&=\prod_{i=1}^n(1-2\lambda_it)^{-1/2}exp\{ \dfrac{\lambda_it}{1-2\lambda_it}(u_i'\mu)^2 \}

\end{align}
$$
Compare $M(t)$ with the moment generating function of $\chi_r(\mu'A\mu)$, we can derive the iff condition.



<u>Proof of Lemma 2</u>

Suppose $\rank(A)=r, \rank(B)=s$. And we decompose the two matrixes by their eigenvalue:
$$
A=\sum_{i=1}^r\lambda_ie_ie_i'\\
B=\sum_{j=1}^s\eta_jf_jf_j'
$$
Then
$$
X'AX=\sum_{i=1}^r\lambda_i(e_i'X)'(e_i'X)\\
X'BX=\sum_{j=1}^s\eta_j(f_j'X)'(f'_jX)
$$
"$\Leftarrow$":

$AB=0$ implies 
$$
\begin{align}
0&=e_i'ABf_j\\
&=\lambda_i\eta_je_i'f_j

\end{align}
$$
which furthermore implies $Cov(e_i'X,f_j'X)=e_i'f_j=0$.

Thus $e_i'X\perp f_j'X$ for all $i,j$.

Finally we conclud $X'AX\perp X'BX$.

"$\Rightarrow$":

For small $t$, 
$$
\begin{align}
K_{X'AX}(t)&=\ln M_{X'AX}(t)\\
&=\sum_{i=1}^n\left\{\dfrac{\lambda_it}{1-2\lambda_it}(e_i'\mu)'(e_i'\mu)-\dfrac12\ln\left( 1-2\lambda_it \right)
\right\}\\
\text{(Taylor Series)}&=\sum_{i=1}^n\dfrac12\sum_{k=1}^{\infty}\left( (e_i'\mu)'(e_i'\mu)\lambda_i^k+\dfrac{\lambda_i^k}{k} \right)(2\lambda_it)^k\\
&=\sum_{k=1}^{\infty}2^{k-1}\left( \mu'A^k\mu+\dfrac{\trace A^k}k \right)t^k
\end{align}
$$
Similarly,
$$
K_{X'BX}(t)=\sum_{k=1}^{\infty}2^{k-1}\left( \mu'B^k\mu+\dfrac{\trace B^k}k \right)t^k
$$
Given $X'AX \perp X'BX$, their moment generating functions are mutipliable.

Thus, for all $\alpha,\beta$,
$$
K_{X'(\alpha A+\beta B)X}=K_{\alpha X'AX}+K_{\beta X'BX}
$$
We can complete the proof by comparing the coeffecients of terms.



<u>Proof of Thm 3</u>

We left this proof in the latter part because it's similar to the former ones.

Without loss of generality, we set $\Sigma=I_n$.

Decompose $B$ by its eigenvalues $B=\sum_{i=1}^n\lambda_iu_iu_i'$.

Thus, $Y'BY=\sum_{i=1}^n\lambda_i(u_i'Y)'(u_i'Y)$.
The given condition $A B=0$ implies $Au_i=0$.
This further implies $Cov(AY, u_i'Y)=0$, which leads to their independence.
The conclusion is trivial then.



Using two theorems we have proved above, we can easily derive some common corollary.

<u>Ex.</u> (Fisher's Lemma)

Suppose $X_i\overset{i.i.d.}{\sim}N(\mu,\sigma^2)$. We then have:

- $\bar X=\dfrac1n\sum_{i=1}^nX_i\sim N(\mu,\dfrac{\sigma^2}{n})$;
- $\dfrac{1}{\sigma^2}\sum_{i=1}^n(\bar X-X_i)^2\sim\chi^2_{n-1}$;
- The above two statistics are independent.

<u>proof</u>

"$1$" is trival.

We now prove "$2$".

Denote $X=(X_1,\dots,X_n)'\sim N(\mu 1_n,\sigma^2I_n)$.

Without loss of generality, we set $\mu=0$.

Notice $\bar X=\dfrac1n1_n'X$, we have
$$
\begin{align}
\dfrac1{\sigma^2}\sum_{i=1}^n(\bar X-X_i)^2&=\dfrac{1}{\sigma^2}(X-\bar X1_n)'(X-\bar X1_n)\\
&=\dfrac{1}{\sigma^2}(X-\dfrac1n1_n1_n'X)'(X-\dfrac1n1_n1_n'X)\\
&=\dfrac{1}{\sigma^2}X'(I_n-\dfrac1n1_n1_n')^2X\\
&=\dfrac{1}{\sigma^2}X'(I_n-\dfrac1n1_n1_n')X

\end{align}
$$
Notice $I_n-\dfrac1n1_n1_n'$ is idempotent with $\rank$ $n-1$.

Thus, we can derive "$2$".

"$3$": The proof is done by checking $(I_n-\dfrac1n1_n1_n')1_n=0$.