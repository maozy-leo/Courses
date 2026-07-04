# Represent Multi-Factor Design Model with Kronecker Product

In multi-factor model, when the design is balanced, the design matrix has a very clean form represented by Kronecker Product.

## Notation and Setup

Consider \(n\) factors \(A_1, A_2, \ldots, A_n\) with \(a_1, a_2, \ldots, a_n\) levels respectively, and \(r\) replications per cell. The total number of observations is

\[
N = a_1 \, a_2 \, \cdots \, a_n \cdot r.
\]

We arrange the response vector \(\mathbf{y}\) so that the last index (replication) varies fastest, and the first factor index varies slowest.

## Building Blocks: Two Column Matrices per Factor

For each factor \(A_i\) with \(a_i\) levels under sum-to-zero constraints, we associate two column matrices:

**The ones vector** \(\mathbf{1}_{a_i}\), an \(a_i \times 1\) vector. It spans the subspace of the grand mean along factor \(i\).

**The contrast matrix** \(C_{a_i}\), an \(a_i \times (a_i - 1)\) matrix satisfying

\[
\mathbf{1}_{a_i}^{\!\top} C_{a_i} = \mathbf{0}^{\!\top}, \qquad \operatorname{rank}(C_{a_i}) = a_i - 1.
\]

The first condition encodes the **sum-to-zero constraint**: every column of \(C_{a_i}\) sums to zero. Together, \(\mathbf{1}_{a_i}\) and \(C_{a_i}\) span all of \(\mathbb{R}^{a_i}\), reflecting the decomposition of the identity into mean and deviation subspaces.

A convenient (though not unique) choice is the set of orthonormal eigenvectors of the centering matrix \(Q_{a_i} = I_{a_i} - \frac{1}{a_i}\mathbf{1}_{a_i}\mathbf{1}_{a_i}^{\!\top}\), so that \(C_{a_i}^{\!\top} C_{a_i} = I_{a_i - 1}\). Other common choices include Helmert contrasts and treatment-vs-baseline contrasts.

## The Additive Model and Its Block Design Matrix

### General \(n\)-factor model

The full factorial effects model under sum-to-zero constraints is

\[
\mathbf{y} = \sum_{S \subseteq \{1,\ldots,n\}} X_S \,\boldsymbol{\theta}_S + \boldsymbol{\epsilon},
\]

where each **column block** \(X_S\) of the design matrix is a Kronecker product:

\[
\boxed{\;
X_S = \bigotimes_{i=1}^{n} Z_i^{(S)} \otimes \mathbf{1}_r, \qquad
Z_i^{(S)} =
\begin{cases}
C_{a_i}, & i \in S, \\[4pt]
\mathbf{1}_{a_i}, & i \notin S.
\end{cases}
\;}
\]

Stacking all \(2^n\) blocks side by side, the full design matrix is

\[
X = \Bigl[\; X_\emptyset \;\Big|\; X_{\{1\}} \;\Big|\; X_{\{2\}} \;\Big|\; \cdots \;\Big|\; X_{\{1,\ldots,n\}} \;\Bigr],
\]

and the full parameter vector is \(\boldsymbol{\theta} = (\boldsymbol{\theta}_\emptyset^{\!\top},\, \boldsymbol{\theta}_{\{1\}}^{\!\top},\, \boldsymbol{\theta}_{\{2\}}^{\!\top},\, \ldots,\, \boldsymbol{\theta}_{\{1,\ldots,n\}}^{\!\top})^{\!\top}\).

### Dimensions of each block

The block \(X_S\) has dimensions \(N \times p_S\), where the number of columns is

\[
p_S = \prod_{i \in S}(a_i - 1),
\]

with the convention that \(p_\emptyset = 1\) (the grand mean is a scalar). One can verify

\[
\sum_{S \subseteq \{1,\ldots,n\}} p_S = \prod_{i=1}^{n}\bigl[1 + (a_i - 1)\bigr] = \prod_{i=1}^{n} a_i,
\]

which equals the total number of cells — confirming that the model is saturated (before counting replication).

## Identification of Effects

Each subset \(S\) indexes a specific factorial effect:

|    Subset \(S\)    | Column Block \(X_S\)                                         | Effect                                 | Parameter                             |                Columns                |
| :----------------: | :----------------------------------------------------------- | :------------------------------------- | :------------------------------------ | :-----------------------------------: |
|   \(\emptyset\)    | \(\mathbf{1}_{a_1} \otimes \cdots \otimes \mathbf{1}_{a_n} \otimes \mathbf{1}_r\) | Grand mean                             | \(\mu\)                               |                 \(1\)                 |
|     \(\{i\}\)      | \(\mathbf{1}_{a_1} \otimes \cdots \otimes C_{a_i} \otimes \cdots \otimes \mathbf{1}_{a_n} \otimes \mathbf{1}_r\) | Main effect of \(A_i\)                 | \(\boldsymbol{\alpha}^{(i)}\)         |              \(a_i - 1\)              |
|    \(\{i,j\}\)     | \(\cdots \otimes C_{a_i} \otimes \cdots \otimes C_{a_j} \otimes \cdots \otimes \mathbf{1}_r\) | Two-way interaction \(A_i \times A_j\) | \(\boldsymbol{\gamma}^{(ij)}\)        |      \((a_i\!-\!1)(a_j\!-\!1)\)       |
|   \(\{i,j,k\}\)    | \(\cdots \otimes C_{a_i} \otimes \cdots \otimes C_{a_j} \otimes \cdots \otimes C_{a_k} \otimes \cdots \otimes \mathbf{1}_r\) | Three-way interaction                  | \(\boldsymbol{\gamma}^{(ijk)}\)       | \((a_i\!-\!1)(a_j\!-\!1)(a_k\!-\!1)\) |
| \(\{1,\ldots,n\}\) | \(C_{a_1} \otimes C_{a_2} \otimes \cdots \otimes C_{a_n} \otimes \mathbf{1}_r\) | Highest-order interaction              | \(\boldsymbol{\gamma}^{(1\cdots n)}\) |        \(\prod_{i}(a_i - 1)\)         |

The pattern is clear: **replacing \(\mathbf{1}_{a_i}\) by \(C_{a_i}\) in position \(i\) "activates" factor \(A_i\) in the effect.** The grand mean activates no factor; a main effect activates exactly one; a \(k\)-way interaction activates exactly \(k\).

## Orthogonality of Blocks

### The key property

**Theorem.** For any two distinct subsets \(S \neq S'\),

\[
X_S^{\!\top}\, X_{S'} = 0.
\]

*Proof.* By the mixed-product rule of Kronecker products,

\[
X_S^{\!\top} X_{S'} = r \cdot \bigotimes_{i=1}^{n} \bigl(Z_i^{(S)}\bigr)^{\!\top} Z_i^{(S')}.
\]

Since \(S \neq S'\), there exists an index \(i\) where one of \(Z_i^{(S)}, Z_i^{(S')}\) equals \(C_{a_i}\) and the other equals \(\mathbf{1}_{a_i}\). Their cross product is

\[
C_{a_i}^{\!\top}\,\mathbf{1}_{a_i} = \mathbf{0} \qquad (\text{the sum-to-zero constraint}),
\]

so the Kronecker product contains a zero factor, and the entire block product vanishes. \(\square\)

### The within-block Gram matrix

For a single block,

\[
X_S^{\!\top} X_S = r \cdot \bigotimes_{i=1}^{n} \bigl(Z_i^{(S)}\bigr)^{\!\top} Z_i^{(S)}
= r \cdot \bigotimes_{i \in S} \bigl(C_{a_i}^{\!\top} C_{a_i}\bigr) \;\otimes\; \bigotimes_{i \notin S} \bigl(\mathbf{1}_{a_i}^{\!\top}\mathbf{1}_{a_i}\bigr)
= r \cdot \Bigl(\prod_{i \notin S} a_i\Bigr) \cdot \bigotimes_{i \in S} C_{a_i}^{\!\top} C_{a_i}.
\]

If orthonormal contrasts are chosen so that \(C_{a_i}^{\!\top} C_{a_i} = I_{a_i - 1}\), this simplifies to

\[
X_S^{\!\top} X_S = r \cdot \Bigl(\prod_{i \notin S} a_i\Bigr) \cdot I_{p_S},
\]

a scalar multiple of the identity — meaning columns within each block are also orthogonal.

## Consequences for OLS Estimation

### Block-diagonal normal equations

The full Gram matrix is

\[
X^{\!\top} X = \bigoplus_{S \subseteq \{1,\ldots,n\}} X_S^{\!\top} X_S,
\]

which is **block-diagonal**. This means the normal equations decouple completely across effects.

### Closed-form OLS estimates

The OLS estimate for effect \(S\) depends only on \(X_S\) and \(\mathbf{y}\):

\[
\hat{\boldsymbol{\theta}}_S = (X_S^{\!\top} X_S)^{-1}\, X_S^{\!\top}\, \mathbf{y}.
\]

With orthonormal contrasts, this becomes

\[
\hat{\boldsymbol{\theta}}_S = \frac{1}{r \,\prod_{i \notin S} a_i} \; X_S^{\!\top}\, \mathbf{y}.
\]

Adding or removing any other effect from the model **does not change** \(\hat{\boldsymbol{\theta}}_S\). This is a hallmark of balanced designs.

### Additive sum-of-squares decomposition

The fitted values for effect \(S\) are \(\hat{\mathbf{y}}_S = X_S \hat{\boldsymbol{\theta}}_S\), and the corresponding sum of squares is

\[
\mathrm{SS}_S = \hat{\boldsymbol{\theta}}_S^{\!\top}\, X_S^{\!\top} X_S \,\hat{\boldsymbol{\theta}}_S = \mathbf{y}^{\!\top} X_S (X_S^{\!\top} X_S)^{-1} X_S^{\!\top}\, \mathbf{y}.
\]

The total decomposition is

\[
\boxed{\;
\|\mathbf{y}\|^2 = \sum_{S \subseteq \{1,\ldots,n\}} \mathrm{SS}_S \;+\; \mathrm{SS}_E,
\;}
\]

where \(\mathrm{SS}_E = \mathbf{y}^{\!\top}(I_{\prod a_i} \otimes Q_r)\,\mathbf{y}\) with \(Q_r = I_r - \frac{1}{r}\mathbf{1}_r\mathbf{1}_r^{\!\top}\) and \(\mathrm{df}_E = \prod_{i=1}^n a_i\,(r-1)\). The \(F\)-statistic for testing effect \(S\) is

\[
F_S = \frac{\mathrm{SS}_S / p_S}{\mathrm{SS}_E / \mathrm{df}_E}.
\]

## Illustrative Example: Two-Factor Model

Consider factors \(A\) (with \(a\) levels) and \(B\) (with \(b\) levels) and \(r\) replicates per cell. The model is

\[
y_{ijk} = \mu + \alpha_i + \beta_j + (\alpha\beta)_{ij} + \epsilon_{ijk},
\]

subject to \(\sum_i \alpha_i = 0\), \(\sum_j \beta_j = 0\), \(\sum_i (\alpha\beta)_{ij} = 0\;\forall\,j\), \(\sum_j (\alpha\beta)_{ij} = 0\;\forall\,i\).

The design matrix partitions into four column blocks:

\[
X = \Bigl[\; \underbrace{\mathbf{1}_a \otimes \mathbf{1}_b \otimes \mathbf{1}_r}_{X_\emptyset} \;\Big|\; \underbrace{C_a \otimes \mathbf{1}_b \otimes \mathbf{1}_r}_{X_A} \;\Big|\; \underbrace{\mathbf{1}_a \otimes C_b \otimes \mathbf{1}_r}_{X_B} \;\Big|\; \underbrace{C_a \otimes C_b \otimes \mathbf{1}_r}_{X_{AB}} \;\Bigr].
\]

Orthogonality holds pairwise among all four blocks. For instance,

\[
X_A^{\!\top} X_B = r \cdot (C_a^{\!\top}\mathbf{1}_a) \otimes (\mathbf{1}_b^{\!\top} C_b) = r \cdot \mathbf{0} \otimes \mathbf{0} = 0.
\]

With orthonormal contrasts, the within-block Gram matrices are

\[
X_\emptyset^{\!\top} X_\emptyset = abr, \quad
X_A^{\!\top} X_A = br \cdot I_{a-1}, \quad
X_B^{\!\top} X_B = ar \cdot I_{b-1}, \quad
X_{AB}^{\!\top} X_{AB} = r \cdot I_{(a-1)(b-1)}.
\]

The OLS estimates recover the familiar marginal-mean formulas:

\[
\hat{\mu} = \bar{y}_{\cdot\cdot\cdot},
\]

\[
\hat{\alpha}_i = \bar{y}_{i\cdot\cdot} - \bar{y}_{\cdot\cdot\cdot},
\]

\[
\hat{\beta}_j = \bar{y}_{\cdot j\cdot} - \bar{y}_{\cdot\cdot\cdot},
\]

\[
\widehat{(\alpha\beta)}_{ij} = \bar{y}_{ij\cdot} - \bar{y}_{i\cdot\cdot} - \bar{y}_{\cdot j\cdot} + \bar{y}_{\cdot\cdot\cdot}.
\]

The ANOVA table is:

| Source                     | \(X_S\)                                                    |          SS          |       df       |
| :------------------------- | :--------------------------------------------------------- | :------------------: | :------------: |
| Grand Mean \(\mu\)         | \(\mathbf{1}_a \otimes \mathbf{1}_b \otimes \mathbf{1}_r\) | \(\mathrm{SS}_\mu\)  |     \(1\)      |
| Main Effect \(A\)          | \(C_a \otimes \mathbf{1}_b \otimes \mathbf{1}_r\)          |  \(\mathrm{SS}_A\)   |   \(a - 1\)    |
| Main Effect \(B\)          | \(\mathbf{1}_a \otimes C_b \otimes \mathbf{1}_r\)          |  \(\mathrm{SS}_B\)   |   \(b - 1\)    |
| Interaction \(A \times B\) | \(C_a \otimes C_b \otimes \mathbf{1}_r\)                   | \(\mathrm{SS}_{AB}\) | \((a-1)(b-1)\) |
| Error                      | —                                                          |  \(\mathrm{SS}_E\)   |  \(ab(r-1)\)   |
| **Total**                  | —                                                          | \(\|\mathbf{y}\|^2\) |    \(abr\)     |

## Summary of Key Properties

The column-block Kronecker representation under sum-to-zero constraints yields the following elegant properties for balanced factorial designs:

1. **Explicit block partition.** The design matrix \(X\) splits into \(2^n\) column blocks, each a single Kronecker product of \(\mathbf{1}\)'s and \(C\)'s, directly encoding which factors participate in each effect.

2. **Inter-block orthogonality.** \(X_S^{\!\top} X_{S'} = 0\) for \(S \neq S'\), guaranteed by the sum-to-zero constraint \(C_{a_i}^{\!\top}\mathbf{1}_{a_i} = \mathbf{0}\). This makes \(X^{\!\top}X\) block-diagonal.

3. **Decoupled estimation.** The OLS estimate for every effect is computed from its own block alone — no other effect interferes. Adding or removing terms does not alter the remaining estimates.

4. **Type I = Type II = Type III SS.** Because of inter-block orthogonality, sequential, partial, and marginal sums of squares all coincide. The order of model fitting is irrelevant.

5. **Intra-block orthogonality (with orthonormal contrasts).** Choosing \(C_{a_i}^{\!\top}C_{a_i} = I_{a_i-1}\) further simplifies each \(X_S^{\!\top}X_S\) to a scalar multiple of the identity, so even individual contrast coefficients within each effect are estimated independently.

6. **Transparent degrees of freedom.** The df for effect \(S\) is \(p_S = \prod_{i \in S}(a_i - 1)\), which is simply the number of columns in \(X_S\) — the product of the ranks of the contrast matrices involved.



# Appendix

A brief review of Krnocker product is attached here.

Let

$$
A=(a_{ij})\in \mathbb{F}^{m\times n},\qquad B\in \mathbb{F}^{p\times q},
$$

where $\mathbb{F}=\mathbb{R}$ or $\mathbb{C}$. The **Kronecker product** is

$$
A\otimes B=
\begin{bmatrix}
a_{11}B & \cdots & a_{1n}B\\
\vdots & \ddots & \vdots\\
a_{m1}B & \cdots & a_{mn}B
\end{bmatrix}\in \mathbb{F}^{mp\times nq}.
$$

I will prove the main properties one by one. The cleanest way is to use the block definition directly.

---

## 1. Size

### Claim

If $A$ is $m\times n$ and $B$ is $p\times q$, then

$$
A\otimes B\in \mathbb{F}^{mp\times nq}.
$$

### Proof

The matrix $A\otimes B$ has $m$ block rows and $n$ block columns. Each block $a_{ij}B$ has size $p\times q$. Therefore the total number of rows is $mp$, and the total number of columns is $nq$. ∎

---

## 2. Bilinearity

### Claim

$$
(A_1+A_2)\otimes B=A_1\otimes B+A_2\otimes B,
$$

$$
A\otimes (B_1+B_2)=A\otimes B_1+A\otimes B_2,
$$

$$
(cA)\otimes B=A\otimes (cB)=c(A\otimes B).
$$

### Proof

Write $A_1=(a_{ij}^{(1)})$, $A_2=(a_{ij}^{(2)})$. Then

$$
(A_1+A_2)\otimes B
=\big[(a_{ij}^{(1)}+a_{ij}^{(2)})B\big]_{ij}.
$$

Since scalar multiplication distributes over matrix addition,

$$
(a_{ij}^{(1)}+a_{ij}^{(2)})B
=a_{ij}^{(1)}B+a_{ij}^{(2)}B.
$$

Hence blockwise,

$$
\begin{aligned}
(A_1+A_2)\otimes B
&=\big[a_{ij}^{(1)}B+a_{ij}^{(2)}B\big]_{ij}\\
&=\big[a_{ij}^{(1)}B\big]_{ij}+\big[a_{ij}^{(2)}B\big]_{ij}\\
&=A_1\otimes B+A_2\otimes B.
\end{aligned}
$$

Similarly,

$$
\begin{aligned}
A\otimes (B_1+B_2)
&=[a_{ij}(B_1+B_2)]_{ij}\\
&=[a_{ij}B_1+a_{ij}B_2]_{ij}\\
&=A\otimes B_1+A\otimes B_2.
\end{aligned}
$$

Finally,

$$
(cA)\otimes B=[(ca_{ij})B]_{ij}=[c(a_{ij}B)]_{ij}=c(A\otimes B),
$$

and similarly

$$
A\otimes (cB)=[a_{ij}(cB)]_{ij}=c(A\otimes B).
$$

∎

---

## 3. Mixed-product rule

### Claim

If $A\in\mathbb{F}^{m\times n}$, $C\in\mathbb{F}^{n\times r}$, $B\in\mathbb{F}^{p\times q}$, $D\in\mathbb{F}^{q\times s}$, then

$$
(A\otimes B)(C\otimes D)=(AC)\otimes (BD).
$$

### Proof

Write $A=(a_{ij})$, $C=(c_{jk})$. Then $A\otimes B$ is an $m\times n$ block matrix with $(i,j)$-block $a_{ij}B$, and $C\otimes D$ is an $n\times r$ block matrix with $(j,k)$-block $c_{jk}D$.

The $(i,k)$-block of $(A\otimes B)(C\otimes D)$ is

$$
\sum_{j=1}^n (a_{ij}B)(c_{jk}D).
$$

Since $c_{jk}$ is a scalar,

$$
(a_{ij}B)(c_{jk}D)=a_{ij}c_{jk}(BD).
$$

Therefore

$$
\sum_{j=1}^n (a_{ij}B)(c_{jk}D)
=\left(\sum_{j=1}^n a_{ij}c_{jk}\right)BD.
$$

But $\sum_{j=1}^n a_{ij}c_{jk}$ is exactly the $(i,k)$-entry of $AC$. Hence the $(i,k)$-block is

$$
(AC)_{ik}\,BD.
$$

So the full block matrix is precisely

$$
(AC)\otimes (BD).
$$

∎

---

## 4. Transpose

### Claim

$$
(A\otimes B)^T=A^T\otimes B^T.
$$

### Proof

The $(i,j)$-block of $A\otimes B$ is $a_{ij}B$. When we transpose a block matrix, blocks swap positions and each block is transposed. So the $(j,i)$-block of $(A\otimes B)^T$ is

$$
(a_{ij}B)^T=a_{ij}B^T.
$$

On the other hand, the $(j,i)$-block of $A^T\otimes B^T$ is

$$
(A^T)_{ji}B^T=a_{ij}B^T.
$$

The blocks agree everywhere, so

$$
(A\otimes B)^T=A^T\otimes B^T.
$$

∎

---

## 5. Conjugate transpose

### Claim

Over $\mathbb{C}$,

$$
(A\otimes B)^*=A^*\otimes B^*.
$$

### Proof

Recall $M^*=\overline{M}^{\,T}$. Then

$$
(A\otimes B)^*
=\overline{A\otimes B}^{\,T}.
$$

Since entrywise complex conjugation distributes over scalar multiplication,

$$
\overline{A\otimes B}=\overline{A}\otimes \overline{B}.
$$

Using the transpose property already proved,

$$
\begin{aligned}
(\overline{A}\otimes \overline{B})^T
&=\overline{A}^{\,T}\otimes \overline{B}^{\,T}\\
&=A^*\otimes B^*.
\end{aligned}
$$

∎

---

## 6. Inverse

### Claim

If $A$ and $B$ are invertible, then

$$
(A\otimes B)^{-1}=A^{-1}\otimes B^{-1}.
$$

### Proof

Use the mixed-product rule:

$$
\begin{aligned}
(A\otimes B)(A^{-1}\otimes B^{-1})
&=(AA^{-1})\otimes (BB^{-1})\\
&=I_m\otimes I_p.
\end{aligned}
$$

Now $I_m\otimes I_p$ is the identity matrix of size $mp\times mp$. Indeed, it is a block diagonal matrix with $m$ diagonal blocks $I_p$. So

$$
I_m\otimes I_p=I_{mp}.
$$

Similarly,

$$
(A^{-1}\otimes B^{-1})(A\otimes B)=I_{nq}
$$

when $A,B$ are square of sizes $n,p$ respectively. Therefore $A^{-1}\otimes B^{-1}$ is the inverse of $A\otimes B$. ∎

---

## 7. Determinant

### Claim

If $A\in\mathbb{F}^{n\times n}$, $B\in\mathbb{F}^{m\times m}$, then

$$
\det(A\otimes B)=\det(A)^m\det(B)^n.
$$

### Proof

I will give a standard proof using eigenvalues.

Assume first that $\mathbb{F}=\mathbb{C}$. Let $\lambda_1,\dots,\lambda_n$ be the eigenvalues of $A$, and $\mu_1,\dots,\mu_m$ the eigenvalues of $B$, counted with algebraic multiplicity. A later property shows that the eigenvalues of $A\otimes B$ are exactly

$$
\lambda_i\mu_j,\qquad 1\le i\le n,\ 1\le j\le m.
$$

Therefore

$$
\det(A\otimes B)
=\prod_{i=1}^n\prod_{j=1}^m (\lambda_i\mu_j).
$$

Rearrange the product:

$$
\begin{aligned}
\prod_{i=1}^n\prod_{j=1}^m (\lambda_i\mu_j)
&=\left(\prod_{i=1}^n \lambda_i^m\right)\left(\prod_{j=1}^m \mu_j^n\right)\\
&=\left(\prod_{i=1}^n \lambda_i\right)^m
\left(\prod_{j=1}^m \mu_j\right)^n.
\end{aligned}
$$

But the determinant equals the product of eigenvalues, so

$$
\det(A\otimes B)=\det(A)^m\det(B)^n.
$$

For real matrices, the same identity follows because the determinant is unchanged whether we view the matrices over $\mathbb{R}$ or $\mathbb{C}$. ∎

### Remark

A proof avoiding eigenvalues is possible but longer. The eigenvalue argument is the cleanest once the eigenvalue property is established.

---

## 8. Trace

### Claim

If $A,B$ are square, then

$$
\operatorname{tr}(A\otimes B)=\operatorname{tr}(A)\operatorname{tr}(B).
$$

### Proof

Let $A=(a_{ij})\in\mathbb{F}^{n\times n}$. The block matrix $A\otimes B$ has diagonal blocks

$$
a_{11}B,\ a_{22}B,\ \dots,\ a_{nn}B.
$$

The trace of a block matrix is the sum of the traces of its diagonal blocks. Hence

$$
\operatorname{tr}(A\otimes B)
=\sum_{i=1}^n \operatorname{tr}(a_{ii}B).
$$

Since $a_{ii}$ is a scalar,

$$
\operatorname{tr}(a_{ii}B)=a_{ii}\operatorname{tr}(B).
$$

So

$$
\begin{aligned}
\operatorname{tr}(A\otimes B)
&=\sum_{i=1}^n a_{ii}\operatorname{tr}(B)\\
&=\left(\sum_{i=1}^n a_{ii}\right)\operatorname{tr}(B)\\
&=\operatorname{tr}(A)\operatorname{tr}(B).
\end{aligned}
$$

∎

---

## 9. Rank

### Claim

$$
\operatorname{rank}(A\otimes B)=\operatorname{rank}(A)\operatorname{rank}(B).
$$

### Proof

This is most transparent using rank normal forms.

Let

$$
r=\operatorname{rank}(A),\qquad s=\operatorname{rank}(B).
$$

There exist invertible matrices $P,Q,R,S$ such that

$$
PAQ=
\begin{bmatrix}
I_r & 0\\
0 & 0
\end{bmatrix},
\qquad
RBS=
\begin{bmatrix}
I_s & 0\\
0 & 0
\end{bmatrix}.
$$

Now use the mixed-product rule:

$$
(P\otimes R)(A\otimes B)(Q\otimes S)
=(PAQ)\otimes (RBS).
$$

Since $P\otimes R$ and $Q\otimes S$ are invertible, multiplication by them does not change rank. Thus

$$
\operatorname{rank}(A\otimes B)
=\operatorname{rank}\big((PAQ)\otimes (RBS)\big).
$$

So it remains to compute the rank of

$$
\begin{bmatrix}
I_r & 0\\
0 & 0
\end{bmatrix}
\otimes
\begin{bmatrix}
I_s & 0\\
0 & 0
\end{bmatrix}.
$$

This is a block diagonal-type matrix whose nonzero part is exactly $I_r\otimes I_s$, and

$$
I_r\otimes I_s=I_{rs}.
$$

Hence its rank is $rs$. Therefore

$$
\operatorname{rank}(A\otimes B)=rs=\operatorname{rank}(A)\operatorname{rank}(B).
$$

∎

---

## 10. Eigenvalues

### Claim

If $x\neq 0$ is an eigenvector of $A$ with eigenvalue $\lambda$, and $y\neq 0$ is an eigenvector of $B$ with eigenvalue $\mu$, then $x\otimes y\neq 0$ is an eigenvector of $A\otimes B$ with eigenvalue $\lambda\mu$. Hence the eigenvalues of $A\otimes B$ are all pairwise products $\lambda_i\mu_j$.

### Proof

First prove the eigenvector statement:

$$
(A\otimes B)(x\otimes y)=(Ax)\otimes (By).
$$

This is a special case of the mixed-product rule, viewing vectors as matrices with one column. Since $Ax=\lambda x$ and $By=\mu y$,

$$
\begin{aligned}
(A\otimes B)(x\otimes y)
&=(\lambda x)\otimes (\mu y)\\
&=\lambda\mu (x\otimes y).
\end{aligned}
$$

Thus $x\otimes y$ is an eigenvector with eigenvalue $\lambda\mu$.

So every product $\lambda_i\mu_j$ is an eigenvalue.

To show these are all eigenvalues, suppose $A$ has size $n\times n$ and $B$ has size $m\times m$. Then $A\otimes B$ has size $nm\times nm$, so it has $nm$ eigenvalues over $\mathbb{C}$, counted with algebraic multiplicity.

If $A$ and $B$ are diagonalizable, choose bases of eigenvectors and write

$$
A=P\Lambda P^{-1},\qquad B=QMQ^{-1},
$$

with $\Lambda=\operatorname{diag}(\lambda_1,\dots,\lambda_n)$, $M=\operatorname{diag}(\mu_1,\dots,\mu_m)$. Then

$$
A\otimes B
=(P\otimes Q)(\Lambda\otimes M)(P^{-1}\otimes Q^{-1}).
$$

So $A\otimes B$ is similar to $\Lambda\otimes M$, which is diagonal with diagonal entries $\lambda_i\mu_j$. Hence those are exactly the eigenvalues.

For general matrices, approximate by diagonalizable matrices or use Schur form; the conclusion remains true with algebraic multiplicity. ∎

### Remark

The key concrete identity is

$$
(A\otimes B)(x\otimes y)=(Ax)\otimes(By).
$$

---

## 11. Positive semidefiniteness / positive definiteness

### Claim

If $A\succeq 0$ and $B\succeq 0$, then $A\otimes B\succeq 0$. If $A\succ 0$ and $B\succ 0$, then $A\otimes B\succ 0$.

### Proof

First assume $A,B$ are Hermitian positive semidefinite. Then by the spectral theorem,

$$
A=\sum_i \lambda_i u_i u_i^*,\qquad B=\sum_j \mu_j v_j v_j^*,
$$

with $\lambda_i,\mu_j\ge 0$.

Now

$$
A\otimes B
=\sum_{i,j}\lambda_i\mu_j\,(u_i u_i^*)\otimes (v_j v_j^*).
$$

Using

$$
(u_i u_i^*)\otimes (v_j v_j^*)=(u_i\otimes v_j)(u_i\otimes v_j)^*,
$$

we get

$$
A\otimes B
=\sum_{i,j}\lambda_i\mu_j\,(u_i\otimes v_j)(u_i\otimes v_j)^*.
$$

Each term is positive semidefinite, and coefficients $\lambda_i\mu_j$ are nonnegative, so $A\otimes B\succeq 0$.

If $A\succ 0$ and $B\succ 0$, then all $\lambda_i>0$ and $\mu_j>0$. Hence every eigenvalue $\lambda_i\mu_j$ of $A\otimes B$ is strictly positive, so $A\otimes B\succ 0$. ∎

---

# The vectorization identity

This was one of the most useful formulas:

$$
\operatorname{vec}(AXB)=(B^T\otimes A)\operatorname{vec}(X).
$$

I will prove it carefully.

Let $X\in\mathbb{F}^{n\times p}$, $A\in\mathbb{F}^{m\times n}$, $B\in\mathbb{F}^{p\times q}$. Then $AXB\in\mathbb{F}^{m\times q}$. Recall that $\operatorname{vec}(X)$ stacks the columns of $X$ one underneath another.

## Step 1: prove $\operatorname{vec}(uv^T)=v\otimes u$

Let $u\in\mathbb{F}^m$, $v\in\mathbb{F}^q$. Then

$$
uv^T=[v_1u\ \ v_2u\ \ \cdots\ \ v_qu].
$$

Therefore

$$
\operatorname{vec}(uv^T)
=\begin{bmatrix}
v_1u\\
v_2u\\
\vdots\\
v_qu
\end{bmatrix}
=v\otimes u.
$$

## Step 2: write $X$ by columns

Let the columns of $X$ be $x_1,\dots,x_p$. Then

$$
XB=
\bigg[\sum_{j=1}^p b_{j1}x_j,\ \dots,\ \sum_{j=1}^p b_{jq}x_j\bigg].
$$

Hence

$$
AXB=
\bigg[\sum_{j=1}^p b_{j1}Ax_j,\ \dots,\ \sum_{j=1}^p b_{jq}Ax_j\bigg].
$$

Vectorizing,

$$
\operatorname{vec}(AXB)
=\begin{bmatrix}
\sum_{j=1}^p b_{j1}Ax_j\\
\sum_{j=1}^p b_{j2}Ax_j\\
\vdots\\
\sum_{j=1}^p b_{jq}Ax_j
\end{bmatrix}.
$$

## Step 3: compute $(B^T\otimes A)\operatorname{vec}(X)$

We have

$$
\operatorname{vec}(X)=
\begin{bmatrix}
x_1\\ x_2\\ \vdots \\ x_p
\end{bmatrix}.
$$

Now $B^T\otimes A$ is the block matrix whose $(k,j)$-block is $(B^T)_{kj}A=b_{jk}A$. So multiplying by $\operatorname{vec}(X)$,

$$
(B^T\otimes A)\operatorname{vec}(X)
=\begin{bmatrix}
\sum_{j=1}^p b_{j1}Ax_j\\
\sum_{j=1}^p b_{j2}Ax_j\\
\vdots\\
\sum_{j=1}^p b_{jq}Ax_j
\end{bmatrix}.
$$

This is exactly $\operatorname{vec}(AXB)$. Therefore

$$
\operatorname{vec}(AXB)=(B^T\otimes A)\operatorname{vec}(X).
$$

∎

---

# Two auxiliary identities that are often used in proofs

These were used implicitly above.

## Identity A

$$
(\alpha B)(\beta D)=(\alpha\beta)(BD)
$$

for scalars $\alpha,\beta$. This is immediate from ordinary matrix multiplication.

## Identity B

$$
(A\otimes B)(x\otimes y)=(Ax)\otimes (By).
$$

### Proof

View $x$ and $y$ as column matrices. Then this is the mixed-product rule with $C=x$, $D=y$:

$$
(A\otimes B)(x\otimes y)=(Ax)\otimes(By).
$$

∎

---

# Important caution: $A\otimes B\neq B\otimes A$ in general

This is not really a "property to prove," but it is worth clarifying. They usually have different sizes unless $A$ and $B$ are square of matching dimensions. Even when both are defined and square, the block layout is different.

However, they are **permutation similar**: there exists a permutation matrix $P$ such that

$$
P^T(A\otimes B)P=B\otimes A.
$$

That is why they have the same eigenvalues.



