<u>Lemma</u>

Suppose $A_i,i=1,\dots,k$ is symmetric, $A_1+\dots+A_k=I_n$. The following statements are equivalent.

- $A_i^2=A_i,i=1,\dots,k$;
- $\\rank(A_1)+\dots+\rank(A_k)=n$;
- $A_iA_j=0,\forall i\ne j$.

<u>Proof</u>

"$1\Rightarrow2$":

$A_i$ is idempotent, then $\rank(A_i)=\trace A_i,\forall i$.

Taking trace of $A_1+\cdots+A_k=I_n$ completes the proof.

"$2\Rightarrow 3$":

Decompose $A_i$ by its eigenvalue $A_i=U_i\Lambda_iU_i',\forall i$, where $U_i$ is an $r_i\times n$ matrix and diagnal elements of $\Lambda_i$ is non-zero.

Let $U=[U_1,\dots,U_k], \Lambda=diag[\Lambda_1,\dots,\Lambda_k]$, which are $n\times n$ matrices.

Then $U\Lambda U'=\sum_{i=1}^kA_i=I_n$, which implies the invertibility of $U$.

Thus $\Lambda U'U=I_n$.

See $\Lambda U'U$ blockwise:

We have $\Lambda_iU_iU_j'=0,\forall i\ne j$, which implies $A_iA_j=0,\forall i\ne j$.

"$3\Rightarrow1$":

$A_i=A_iI_n=A_i(A_1+\dots+A_k)=A_i^2,\forall i$.



<u>Thm</u> *(Cochran's Theorem)*

Suppose $X\sim N(\mu,I_n)$. $Q_i=X'A_iX$, satisfying $Q_1+\dots+Q_k=X' X$. $A_i$ symmetric and $\rank(A_i)=r_i,i=1,\dots,k$.  Then the following statements are equivalent:

- $Q_i\sim\chi^2_{m_i}(\lambda_i)$, where $\lambda_i=\mu'A_i\mu$;
- $r_1+\dots+r_k=n$;
- $Q_i\perp Q_j,\forall i\ne j$.

<u>proof</u>

$Q_1+\dots+Q_k=X'(A_1+\dots+A_k)X=X'X\sim \chi^2_n(\mu'\mu)$.

$\Rightarrow$ $\sum_{i=1}^nA_i$ is idempotent with $n$ (using chi-square lemma), impling $\sum_{i=1}^nA_i=I_n$

$1\Rightarrow2$:

$Q_i\sim \chi_{m_i}^2(\lambda_i)$ $\Rightarrow$ $A_i$ is idempotent with $\rank$ $m_i$.

$\Rightarrow$$r_i=m_i$

$\Rightarrow$$r_1+\dots+r_k=\trace A_1\dots+\trace A_k=n$.

"$2\Rightarrow3$":

It suffices to show that $A_i A_j=0,\forall i\ne j$, which comes directly from our lemma.

"$3\Rightarrow1$":

Since $Q_i\perp Q_j$, then $A_iA_j=0$.

$\Rightarrow A_i$ is idempotent.

$\Rightarrow$ $Q_i\sim\chi^2_{r_i}(\lambda_i)$.



Notice: We can easily generalize the theorem to the case $X\sim N(\mu,\Sigma)$, using the transformation $Y=\Sigma^{-1/2}X$.
