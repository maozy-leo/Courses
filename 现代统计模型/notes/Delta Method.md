# Delta Method

## Preliminaries

- $X_n \in \mathbb{R}^k$ is a random vector.
- $\theta \in \mathbb{R}^k$ is a constant vector.
- $a_n \to \infty$ is a sequence of positive numbers.
- $g : \mathbb{R}^k \to \mathbb{R}^m$ is a function that is sufficiently smooth near $\theta$.

We will repeatedly use the following notation:

- $Dg(\theta)$ denotes the derivative (Jacobian) of $g$ at $\theta$.
- If $g : \mathbb{R}^k \to \mathbb{R}$, then $Dg(\theta)$ is a $1 \times k$ row vector, equivalently $\nabla g(\theta)'$.
- If $g : \mathbb{R}^k \to \mathbb{R}$ is twice differentiable, then $H_g(\theta)$ denotes its Hessian matrix at $\theta$.

We also use the standard asymptotic notation:

- $o_p(1)$ means convergence to $0$ in probability.
- $O_p(1)$ means bounded in probability.
- $\xrightarrow{d}$ means convergence in distribution.
- $\xrightarrow{p}$ means convergence in probability.

---

## Key tools

### Slutsky's theorem

If $Y_n \xrightarrow{d} Y$ and $Z_n \xrightarrow{p} c$, then

$$
Y_n + Z_n \xrightarrow{d} Y + c
$$

and

$$
Z_n Y_n \xrightarrow{d} cY.
$$

### Continuous mapping theorem

If $Y_n \xrightarrow{d} Y$ and $f$ is continuous, then

$$
f(Y_n) \xrightarrow{d} f(Y).
$$

### Differentiability expansions

If $g$ is differentiable at $\theta$, then

$$
g(\theta + h) = g(\theta) + Dg(\theta)h + r(h),
$$

where

$$
\frac{\|r(h)\|}{\|h\|} \to 0 \quad \text{as } h \to 0.
$$

If $g$ is twice differentiable at $\theta$ and scalar-valued, then

$$
g(\theta + h) = g(\theta) + Dg(\theta)h + \frac{1}{2} h' H_g(\theta) h + r_2(h),
$$

where

$$
\frac{r_2(h)}{\|h\|^2} \to 0 \quad \text{as } h \to 0.
$$

## First-order delta method

### Scalar-valued version

Let $g : \mathbb{R}^k \to \mathbb{R}$. Assume

$$
a_n (X_n - \theta) \xrightarrow{d} Z,
$$

where $Z$ is an $\mathbb{R}^k$-valued random vector. If $g$ is differentiable at $\theta$, then

$$
a_n \bigl(g(X_n) - g(\theta)\bigr) \xrightarrow{d} Dg(\theta) Z = \nabla g(\theta)' Z.
$$

In particular, if

$$
\sqrt{n}(X_n - \theta) \xrightarrow{d} N_k(0, \Sigma),
$$

then

$$
\sqrt{n}\bigl(g(X_n) - g(\theta)\bigr)
\xrightarrow{d}
N\bigl(0, \nabla g(\theta)' \Sigma \nabla g(\theta)\bigr).
$$

### Proof

By differentiability of $g$ at $\theta$,

$$
g(X_n) - g(\theta) = Dg(\theta)(X_n - \theta) + r_n,
$$

with

$$
\frac{r_n}{\|X_n - \theta\|} \xrightarrow{p} 0.
$$

Multiplying by $a_n$ gives

$$
a_n\bigl(g(X_n) - g(\theta)\bigr)
=
Dg(\theta)\, a_n(X_n - \theta) + a_n r_n.
$$

So it remains to show that

$$
a_n r_n \xrightarrow{p} 0.
$$

Write

$$
r_n = \|X_n - \theta\| \eta_n,
\qquad
\eta_n \xrightarrow{p} 0.
$$

Then

$$
a_n r_n = \bigl(a_n \|X_n - \theta\|\bigr)\eta_n.
$$

Since

$$
a_n(X_n - \theta) \xrightarrow{d} Z,
$$

we have

$$
a_n(X_n - \theta) = O_p(1),
$$

hence

$$
a_n \|X_n - \theta\| = O_p(1).
$$

Therefore,

$$
a_n r_n = O_p(1)\, o_p(1) = o_p(1).
$$

So

$$
a_n\bigl(g(X_n) - g(\theta)\bigr)
=
Dg(\theta)\, a_n(X_n - \theta) + o_p(1).
$$

By Slutsky's theorem,

$$
a_n\bigl(g(X_n) - g(\theta)\bigr) \xrightarrow{d} Dg(\theta) Z.
$$

## First-order delta method: vector-valued version

Let $g : \mathbb{R}^k \to \mathbb{R}^m$. Assume

$$
a_n (X_n - \theta) \xrightarrow{d} Z.
$$

If $g$ is differentiable at $\theta$, then

$$
a_n\bigl(g(X_n) - g(\theta)\bigr) \xrightarrow{d} Dg(\theta) Z.
$$

If moreover

$$
\sqrt{n}(X_n - \theta) \xrightarrow{d} N_k(0, \Sigma),
$$

then

$$
\sqrt{n}\bigl(g(X_n) - g(\theta)\bigr)
\xrightarrow{d}
N_m\bigl(0, Dg(\theta)\Sigma Dg(\theta)'\bigr).
$$

### Proof sketch

The argument is identical to the scalar-valued case. We write

$$
g(X_n) - g(\theta) = Dg(\theta)(X_n - \theta) + r_n,
$$

with

$$
\frac{\|r_n\|}{\|X_n - \theta\|} \xrightarrow{p} 0.
$$

Then

$$
a_n\bigl(g(X_n) - g(\theta)\bigr)
=
Dg(\theta)\, a_n(X_n - \theta) + a_n r_n,
$$

and again $a_n r_n = o_p(1)$.

## Equivalent asymptotic linearization

A common equivalent formulation of the first-order delta method is

$$
a_n\bigl(g(X_n) - g(\theta)\bigr) - Dg(\theta)\, a_n(X_n - \theta) \xrightarrow{p} 0.
$$

That is,

$$
a_n\bigl(g(X_n) - g(\theta)\bigr)
=
Dg(\theta)\, a_n(X_n - \theta) + o_p(1).
$$

This shows that, at first order, the transformation $g$ is asymptotically linear.

## Second-order delta method

### Scalar-valued version with vanishing first derivative

Let $g : \mathbb{R}^k \to \mathbb{R}$. Assume

$$
a_n(X_n - \theta) \xrightarrow{d} Z,
$$

and suppose that $g$ is twice continuously differentiable near $\theta$, with

$$
Dg(\theta) = 0.
$$

Then

$$
a_n^2 \bigl(g(X_n) - g(\theta)\bigr)
\xrightarrow{d}
\frac{1}{2} Z' H_g(\theta) Z.
$$

In particular, if

$$
\sqrt{n}(X_n - \theta) \xrightarrow{d} N_k(0, \Sigma),
$$

then

$$
n \bigl(g(X_n) - g(\theta)\bigr)
\xrightarrow{d}
\frac{1}{2} Z' H_g(\theta) Z,
\qquad
Z \sim N_k(0, \Sigma).
$$

### Proof

By the second-order Taylor expansion,

$$
g(X_n) - g(\theta)
=
Dg(\theta)(X_n - \theta)
+
\frac{1}{2}(X_n - \theta)' H_g(\theta)(X_n - \theta)
+
r_n,
$$

where

$$
\frac{r_n}{\|X_n - \theta\|^2} \xrightarrow{p} 0.
$$

Since $Dg(\theta) = 0$, this reduces to

$$
g(X_n) - g(\theta)
=
\frac{1}{2}(X_n - \theta)' H_g(\theta)(X_n - \theta)
+
r_n.
$$

Multiplying by $a_n^2$ gives

$$
a_n^2 \bigl(g(X_n) - g(\theta)\bigr)
=
\frac{1}{2}
\bigl(a_n(X_n - \theta)\bigr)'
H_g(\theta)
\bigl(a_n(X_n - \theta)\bigr)
+
a_n^2 r_n.
$$

Now write

$$
r_n = \|X_n - \theta\|^2 \eta_n,
\qquad
\eta_n \xrightarrow{p} 0.
$$

Then

$$
a_n^2 r_n
=
\bigl(a_n \|X_n - \theta\|\bigr)^2 \eta_n.
$$

Because $a_n(X_n - \theta) = O_p(1)$, we have

$$
\bigl(a_n \|X_n - \theta\|\bigr)^2 = O_p(1),
$$

and hence

$$
a_n^2 r_n = o_p(1).
$$

Therefore,

$$
a_n^2 \bigl(g(X_n) - g(\theta)\bigr)
=
\frac{1}{2}
\bigl(a_n(X_n - \theta)\bigr)'
H_g(\theta)
\bigl(a_n(X_n - \theta)\bigr)
+
o_p(1).
$$

Since the map

$$
x \mapsto \frac{1}{2} x' H_g(\theta)x
$$

is continuous, the continuous mapping theorem implies

$$
\frac{1}{2}
\bigl(a_n(X_n - \theta)\bigr)'
H_g(\theta)
\bigl(a_n(X_n - \theta)\bigr)
\xrightarrow{d}
\frac{1}{2} Z' H_g(\theta) Z.
$$

Finally, Slutsky's theorem gives the result.

## Second-order expansion without assuming $Dg(\theta)=0$

If $g$ is twice continuously differentiable near $\theta$, then

$$
g(X_n) - g(\theta)
=
Dg(\theta)(X_n - \theta)
+
\frac{1}{2}(X_n - \theta)' H_g(\theta)(X_n - \theta)
+
o_p(\|X_n - \theta\|^2).
$$

Equivalently,

$$
a_n^2 \Bigl(g(X_n) - g(\theta) - Dg(\theta)(X_n - \theta)\Bigr)
\xrightarrow{d}
\frac{1}{2} Z' H_g(\theta) Z.
$$

This shows that:

- the first-order term appears at scale $a_n$;
- after removing the linear term, the quadratic term appears at scale $a_n^2$.

## Vector-valued second-order delta method

Let $g : \mathbb{R}^k \to \mathbb{R}^m$, and write

$$
g = (g_1, \dots, g_m)'.
$$

Assume that each component $g_j$ is twice continuously differentiable near $\theta$, and that

$$
a_n(X_n - \theta) \xrightarrow{d} Z.
$$

Then

$$
a_n^2 \Bigl(g(X_n) - g(\theta) - Dg(\theta)(X_n - \theta)\Bigr)
\xrightarrow{d}
\frac{1}{2}
\begin{pmatrix}
Z' H_{g_1}(\theta) Z \\
\vdots \\
Z' H_{g_m}(\theta) Z
\end{pmatrix}.
$$

If in addition $Dg(\theta) = 0$, then

$$
a_n^2 \bigl(g(X_n) - g(\theta)\bigr)
\xrightarrow{d}
\frac{1}{2}
\begin{pmatrix}
Z' H_{g_1}(\theta) Z \\
\vdots \\
Z' H_{g_m}(\theta) Z
\end{pmatrix}.
$$

### Proof sketch

Apply the scalar second-order Taylor expansion to each component $g_j$:

$$
g_j(X_n) - g_j(\theta) - Dg_j(\theta)(X_n - \theta)
=
\frac{1}{2}(X_n - \theta)' H_{g_j}(\theta)(X_n - \theta)
+
r_{n,j},
$$

with

$$
r_{n,j} = o_p(\|X_n - \theta\|^2).
$$

Then stack the components together.