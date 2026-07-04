# Modes of Convergence of Random Variables

## 1. Setup and notation

Let $(\Omega,\mathcal F,\mathbb P)$ be a probability space. Unless stated otherwise, all random variables are real-valued.

We write

- $X_n \xrightarrow{\text{a.s.}} X$ for almost sure convergence,
- $X_n \xrightarrow{\mathbb P} X$ for convergence in probability,
- $X_n \xrightarrow{L^p} X$ for convergence in $L^p$,
- $X_n \xrightarrow{d} X$ for convergence in distribution.

We also write $\mathcal L(X)$ for the law (distribution) of $X$.

A random variable $Y$ is a **version** of $X$ if $X=Y$ almost surely, that is,
$$
\mathbb P(X=Y)=1.
$$
All of the convergence notions below are unchanged if we replace a random variable by a version of it.

A basic but important distinction is this:

- **Almost sure convergence**, **convergence in probability**, and **$L^p$ convergence** compare $X_n$ and $X$ on the same probability space, because they involve the difference $X_n-X$.
- **Convergence in distribution** depends only on the laws $\mathcal L(X_n)$ and $\mathcal L(X)$, so the variables do **not** need to live on a common probability space.

This is one reason weak convergence is especially convenient in asymptotic statistics.

---

## 2. Almost sure convergence

### Definition

We say that $X_n$ converges **almost surely** to $X$ if
$$
\mathbb P\big(\{\omega\in\Omega : X_n(\omega)\to X(\omega)\}\big)=1.
$$
Equivalently, for almost every outcome $\omega$, the numerical sequence $X_n(\omega)$ converges to $X(\omega)$.

### Equivalent formulations

For $\varepsilon>0$, write
$$
A_n(\varepsilon)=\{|X_n-X|>\varepsilon\}.
$$
Then
$$
X_n\xrightarrow{\text{a.s.}}X
$$
if and only if for every $\varepsilon>0$,
$$
\mathbb P\big(A_n(\varepsilon)\ \text{i.o.}\big)=0,
$$
where “i.o.” means “infinitely often,” that is,
$$
A_n(\varepsilon)\ \text{i.o.}
=
\limsup_{n\to\infty} A_n(\varepsilon)
=
\bigcap_{m=1}^\infty \bigcup_{n\ge m} A_n(\varepsilon).
$$

Another very useful identity is
$$
\{X_n \not\to X\}
=
\bigcup_{m=1}^\infty \limsup_{n\to\infty}\{|X_n-X|>1/m\}.
$$
Hence,
$$
X_n\xrightarrow{\text{a.s.}}X
\quad\Longleftrightarrow\quad
\forall \varepsilon>0,\ \mathbb P\big(\limsup_{n} \{|X_n-X|>\varepsilon\}\big)=0.
$$

Use Borel-Cantelli lemma, the above condition can be interpreted as the following theorem:

If for every $\varepsilon>0$, the series is convergent
$$
\sum_{n=1}^{\infty}\mathbb P(|X_n-X|>\varepsilon)<\infty
$$
then, $X_n\overset{a.s.}\to X$.

### Interpretation

Almost sure convergence is the strongest of the four classical modes in the pathwise sense: along almost every sample point, the sequence eventually behaves like an ordinary convergent sequence of real numbers.

### Example

On $([0,1],\mathcal B,\lambda)$, let
$$
X_n(\omega)=\mathbf 1_{(0,1/n)}(\omega).
$$
Then for every $\omega>0$, $X_n(\omega)=0$ for all sufficiently large $n$, so
$$
X_n(\omega)\to 0
$$
for every $\omega\in(0,1]$, and therefore
$$
X_n\xrightarrow{\text{a.s.}}0.
$$

### Example showing a.s. does not imply $L^1$

Again on $([0,1],\mathcal B,\lambda)$, define
$$
X_n(\omega)=n\,\mathbf 1_{(0,1/n)}(\omega).
$$
Then $X_n(\omega)\to 0$ for every $\omega>0$, so $X_n\to 0$ almost surely. But
$$
\mathbb E|X_n|
=
n\cdot \lambda(0,1/n)
=
1,
$$
so $X_n\not\to 0$ in $L^1$.

This example is one of the standard counterexamples showing that almost sure convergence and $L^1$ convergence are not comparable in general.

---

## 3. Convergence in probability

### Definition

We say that $X_n$ converges **in probability** to $X$ if for every $\varepsilon>0$,
$$
\mathbb P(|X_n-X|>\varepsilon)\to 0.
$$

### Equivalent formulations

The following are equivalent:

1. $X_n\xrightarrow{\mathbb P}X$.
2. For every $\varepsilon>0$,
   $$
   \mathbb P(|X_n-X|\ge \varepsilon)\to 0.
   $$
3. $X_n\to X$ in measure on the probability space $(\Omega,\mathcal F,\mathbb P)$.
4. The metric
   $$
   d(U,V)=\mathbb E[\min\{|U-V|,1\}]
   $$
   satisfies $d(X_n,X)\to 0$.

The equivalence of (1) and (3) is just the finite-measure-space fact that “convergence in measure” becomes “convergence in probability” when the total measure is $1$.

### Subsequence characterization

A fundamental theorem is:

**Theorem (subsequence principle).**
$$
X_n\xrightarrow{\mathbb P}X
\quad\Longleftrightarrow\quad
\text{every subsequence }(X_{n_k})\text{ has a further subsequence }(X_{n_{k_j}})
\text{ such that }X_{n_{k_j}}\xrightarrow{\text{a.s.}}X.
$$

**Proof.**

- Suppose $X_n\xrightarrow{\mathbb P}X$, and let $(X_{n_k})$ be a subsequence. Choose a further subsequence $(X_{n_{k_j}})$ so that
  $$
  \mathbb P\big(|X_{n_{k_j}}-X|>2^{-j}\big)<2^{-j}.
  $$
  Then
  $$
  \sum_{j=1}^\infty \mathbb P\big(|X_{n_{k_j}}-X|>2^{-j}\big)<\infty.
  $$
  By the first Borel–Cantelli lemma,
  $$
  |X_{n_{k_j}}-X|>2^{-j}\ \text{i.o.}
  $$
  has probability $0$, which implies $X_{n_{k_j}}\to X$ almost surely.

- Conversely, suppose every subsequence has a further almost surely convergent subsubsequence, but $X_n\not\xrightarrow{\mathbb P}X$. Then there exist $\varepsilon,\delta>0$ and a subsequence $(X_{n_k})$ such that
  $$
  \mathbb P(|X_{n_k}-X|>\varepsilon)>\delta
  $$
  for all $k$. No further subsubsequence can converge almost surely to $X$, because almost sure convergence implies convergence in probability. Contradiction. $\square$

### Example: in probability but not almost surely

Let $X_1,X_2,\dots$ be independent Bernoulli random variables with
$$
\mathbb P(X_n=1)=\frac1n,\qquad \mathbb P(X_n=0)=1-\frac1n.
$$
Then for any $\varepsilon\in(0,1)$,
$$
\mathbb P(|X_n|>\varepsilon)=\mathbb P(X_n=1)=\frac1n\to 0,
$$
so $X_n\xrightarrow{\mathbb P}0$.

But
$$
\sum_{n=1}^\infty \mathbb P(X_n=1)=\sum_{n=1}^\infty \frac1n=\infty,
$$
and by the second Borel–Cantelli lemma (using independence),
$$
\mathbb P(X_n=1\ \text{i.o.})=1.
$$
Hence $X_n$ does **not** converge to $0$ almost surely.

---

## 4. Convergence in \(L^p\)

### Definition

Let $0<p<\infty$. We say that $X_n$ converges to $X$ **in $L^p$** if
$$
\mathbb E|X_n-X|^p\to 0.
$$
For $p=1$ this is often called **convergence in mean**, and for $p=2$ it is called **mean-square convergence**.

For $p\ge 1$, this is equivalent to
$$
\|X_n-X\|_p\to 0,
\qquad
\|Y\|_p=(\mathbb E|Y|^p)^{1/p},
$$
so $L^p$ convergence is convergence in norm.

For $0<p<1$, $\|Y\|_p$ is not a norm, but the convergence notion $\mathbb E|X_n-X|^p\to 0$ still makes perfect sense.

### Relation between \(L^q\) and \(L^p\)

On a probability space, if $0<p\le q<\infty$, then
$$
\|Y\|_p\le \|Y\|_q.
$$
Therefore,
$$
X_n\xrightarrow{L^q}X
\quad\Longrightarrow\quad
X_n\xrightarrow{L^p}X.
$$

So stronger moment convergence implies weaker moment convergence.

### Example: \(L^1\) but not \(L^2\)

On $([0,1],\mathcal B,\lambda)$, let
$$
X_n(\omega)=\sqrt n\,\mathbf 1_{(0,1/n)}(\omega).
$$
Then
$$
\mathbb E|X_n|
=
\sqrt n\cdot \frac1n
=
\frac1{\sqrt n}\to 0,
$$
so $X_n\to 0$ in $L^1$.

But
$$
\mathbb E(X_n^2)
=
n\cdot \frac1n
=
1,
$$
so $X_n\not\to 0$ in $L^2$.

### Example: \(L^p\) need not imply almost sure convergence

Return to the independent Bernoulli example with $\mathbb P(X_n=1)=1/n$. For every $p>0$,
$$
\mathbb E|X_n|^p=\frac1n\to 0,
$$
so
$$
X_n\xrightarrow{L^p}0
\qquad \text{for every }p>0.
$$
But, as shown above, $X_n\not\to 0$ almost surely.

Thus $L^p$ convergence and almost sure convergence are not comparable in general.

---

## 5. Convergence in distribution

Convergence in distribution is also called **weak convergence**.

### Definition via distribution functions

Let
$$
F_{X_n}(x)=\mathbb P(X_n\le x),\qquad F_X(x)=\mathbb P(X\le x).
$$
We say that $X_n$ converges **in distribution** to $X$ if
$$
F_{X_n}(x)\to F_X(x)
\quad
\text{for every continuity point }x\text{ of }F_X.
$$
We write
$$
X_n\xrightarrow{d}X.
$$

The continuity-point condition is essential because distribution functions can jump.

### Equivalent formulations on \(\mathbb R\)

The following are equivalent:

1. $X_n\xrightarrow{d}X$.
2. For every bounded contieFor every closed set $F\subset\mathbb R$,
   $$
   \limsup_{n\to\infty}\mathbb P(X_n\in F)\le \mathbb P(X\in F).
   $$
4. For every open set $G\subset\mathbb R$,
   $$
   \liminf_{n\to\infty}\mathbb P(X_n\in G)\ge \mathbb P(X\in G).
   $$
5. For every Borel set $B$ such that $\mathbb P(X\in\partial B)=0$,
   $$
   \mathbb P(X_n\in B)\to \mathbb P(X\in B).
   $$

These are standard forms of the **Portmanteau theorem**.

### Characteristic functions

The characteristic function of $X$ is
$$
\varphi_X(t)=\mathbb E[e^{itX}],\qquad t\in\mathbb R.
$$
If $X_n\xrightarrow{d}X$, then
$$
\varphi_{X_n}(t)\to \varphi_X(t)
\qquad \text{for every }t\in\mathbb R.
$$

Conversely, by **Lévy’s continuity theorem**, if $\varphi_{X_n}(t)\to \varphi(t)$ for every $t$, and $\varphi$ is continuous at $0$, then $\varphi$ is the characteristic function of some random variable $X$, and
$$
X_n\xrightarrow{d}X.
$$

### Convergence to a constant

A very important special case is:

**Proposition.** For a constant $c$,
$$
X_n\xrightarrow{d}c
\quad\Longleftrightarrow\quad
X_n\xrightarrow{\mathbb P}c.
$$

**Proof.**

- If $X_n\xrightarrow{\mathbb P}c$, then $X_n\xrightarrow{d}c$ because convergence in probability implies convergence in distribution.
- Conversely, if $X_n\xrightarrow{d}c$, then for any $\varepsilon>0$,
  $$
  B_\varepsilon=(c-\varepsilon,c+\varepsilon)
  $$
  is a continuity set for the point mass at $c$, so
  $$
  \mathbb P(|X_n-c|<\varepsilon)=\mathbb P(X_n\in B_\varepsilon)\to 1.
  $$
  Hence $\mathbb P(|X_n-c|>\varepsilon)\to 0$, which is exactly $X_n\xrightarrow{\mathbb P}c$. $\square$

### Example: in distribution but not in probability

Let $X$ be a Rademacher random variable:
$$
\mathbb P(X=1)=\mathbb P(X=-1)=\frac12.
$$
Define
$$
X_n=(-1)^nX.
$$
Then each $X_n$ has the same distribution as $X$, so
$$
X_n\xrightarrow{d}X.
$$
But $X_n$ does not converge in probability to any random variable.

Indeed, the even subsequence is constantly $X$, while the odd subsequence is constantly $-X$. If $X_n\to Y$ in probability, then uniqueness of the probability limit would force $Y=X$ almost surely from the even subsequence and $Y=-X$ almost surely from the odd subsequence, impossible because $\mathbb P(X\neq 0)=1$.

This shows that convergence in distribution is strictly weaker than convergence in probability.

---

## 6. Other useful notions

### Complete convergence

We say that $X_n$ converges **completely** to $X$ if for every $\varepsilon>0$,
$$
\sum_{n=1}^\infty \mathbb P(|X_n-X|>\varepsilon)<\infty.
$$
Then by Borel–Cantelli,
$$
\mathbb P(|X_n-X|>\varepsilon\ \text{i.o.})=0
$$
for each $\varepsilon>0$, so complete convergence implies almost sure convergence.

This is stronger than almost sure convergence. The spike example
$$
X_n=n\,\mathbf 1_{(0,1/n)}
$$
converges almost surely to $0$, but for $\varepsilon=\tfrac12$,
$$
\sum_{n=1}^\infty \mathbb P(|X_n|>\tfrac12)
=
\sum_{n=1}^\infty \frac1n
=
\infty,
$$
so the convergence is not complete.

### Convergence in measure

On a general measure space, convergence in measure means
$$
\mu(|f_n-f|>\varepsilon)\to 0
\qquad \forall \varepsilon>0.
$$
On a probability space this is exactly convergence in probability.

### Almost uniform convergence

We say that $X_n\to X$ **almost uniformly** if for every $\eta>0$ there exists a set $A\in\mathcal F$ with
$$
\mathbb P(A^c)<\eta
$$
such that the convergence is uniform on $A$:
$$
\sup_{\omega\in A}|X_n(\omega)-X(\omega)|\to 0.
$$

Almost uniform convergence implies almost sure convergence. Conversely, on a probability space, **Egorov’s theorem** says that almost sure convergence implies almost uniform convergence. Thus for real-valued random variables on probability spaces, these two notions are essentially equivalent.

---

## 7. Relationship between the different modes of convergence

### Main implication diagram

For $0<p\le q<\infty$,
$$
X_n\xrightarrow{L^q}X
\Longrightarrow
X_n\xrightarrow{L^p}X
\Longrightarrow
X_n\xrightarrow{\mathbb P}X
\Longrightarrow
X_n\xrightarrow{d}X.
$$

Also,
$$
X_n\ \text{converges completely to }X
\Longrightarrow
X_n\xrightarrow{\text{a.s.}}X
\Longrightarrow
X_n\xrightarrow{\mathbb P}X
\Longrightarrow
X_n\xrightarrow{d}X.
$$

On probability spaces,
$$
X_n\xrightarrow{\text{a.s.}}X
\quad\Longleftrightarrow\quad
X_n\to X \text{ almost uniformly}.
$$

Finally,
$$
X_n\xrightarrow{d}c
\quad\Longleftrightarrow\quad
X_n\xrightarrow{\mathbb P}c
$$
for every constant $c$.

### Proofs of the standard implications

#### Almost sure \(\Rightarrow\) probability

Fix $\varepsilon>0$. If $X_n\to X$ almost surely, then
$$
\mathbf 1_{\{|X_n-X|>\varepsilon\}}\to 0
\quad \text{almost surely}.
$$
These indicators are bounded by $1$, so the dominated convergence theorem gives
$$
\mathbb P(|X_n-X|>\varepsilon)
=
\mathbb E\big[\mathbf 1_{\{|X_n-X|>\varepsilon\}}\big]
\to 0.
$$

#### \(L^p \Rightarrow\) probability

For any $\varepsilon>0$, Markov’s inequality yields
$$
\mathbb P(|X_n-X|>\varepsilon)
\le
\frac{\mathbb E|X_n-X|^p}{\varepsilon^p}.
$$
Thus $X_n\xrightarrow{L^p}X$ implies $X_n\xrightarrow{\mathbb P}X$.

#### Probability \(\Rightarrow\) distribution

A clean proof uses the subsequence principle. Suppose $X_n\xrightarrow{\mathbb P}X$. Let $f$ be bounded and continuous. Every subsequence $(X_{n_k})$ has a further almost surely convergent subsubsequence $(X_{n_{k_j}})$ with
$$
X_{n_{k_j}}\xrightarrow{\text{a.s.}}X.
$$
Therefore
$$
f(X_{n_{k_j}})\xrightarrow{\text{a.s.}} f(X),
$$
and since $f$ is bounded, dominated convergence gives
$$
\mathbb E[f(X_{n_{k_j}})]\to \mathbb E[f(X)].
$$
Hence every subsequence of $\mathbb E[f(X_n)]$ has a further subsubsequence converging to $\mathbb E[f(X)]$, which forces
$$
\mathbb E[f(X_n)]\to \mathbb E[f(X)].
$$
By the bounded-continuous-function characterization, $X_n\xrightarrow{d}X$.

#### \(L^q \Rightarrow L^p\) for \(q\ge p\)

On a probability space,
$$
\|Y\|_p\le \|Y\|_q.
$$
Applying this to $Y=X_n-X$ gives
$$
\|X_n-X\|_p\le \|X_n-X\|_q.
$$
So if $\|X_n-X\|_q\to 0$, then $\|X_n-X\|_p\to 0$.

#### Complete convergence \(\Rightarrow\) almost sure convergence

If for every $\varepsilon>0$,
$$
\sum_{n=1}^\infty \mathbb P(|X_n-X|>\varepsilon)<\infty,
$$
then by Borel–Cantelli
$$
\mathbb P(|X_n-X|>\varepsilon\ \text{i.o.})=0
$$
for each $\varepsilon>0$. Applying this to $\varepsilon=1/m$, $m\ge 1$, yields $X_n\to X$ almost surely.

### Which converses fail?

Most converses fail. In particular:

- Convergence in probability does **not** imply almost sure convergence.
- Convergence in distribution does **not** imply convergence in probability.
- Almost sure convergence does **not** imply $L^1$ convergence.
- $L^1$ convergence does **not** imply almost sure convergence.
- $L^1$ convergence does **not** imply $L^2$ convergence.

So the main chains of implications are genuinely one-way.

### Counterexamples separating the modes

1. **In probability but not almost surely:**  
   Independent Bernoulli variables with $\mathbb P(X_n=1)=1/n$.

2. **In distribution but not in probability:**  
   $X_n=(-1)^nX$ with $X$ Rademacher.

3. **Almost surely but not in \(L^1\):**  
   $X_n=n\,\mathbf 1_{(0,1/n)}$ on $([0,1],\mathcal B,\lambda)$.

4. **In \(L^1\) but not in \(L^2\):**  
   $X_n=\sqrt n\,\mathbf 1_{(0,1/n)}$.

5. **In \(L^p\) for every \(p>0\), but not almost surely:**  
   Independent Bernoulli variables with $\mathbb P(X_n=1)=1/n$.

### Important special cases where converses do hold

Some converses become true under extra hypotheses.

#### Weak convergence to a constant

As proved above,
$$
X_n\xrightarrow{d}c
\quad\Longleftrightarrow\quad
X_n\xrightarrow{\mathbb P}c.
$$

#### Uniform integrability upgrades probability to \(L^1\)

A family $\{Y_n\}$ is **uniformly integrable** if
$$
\lim_{K\to\infty}\sup_n \mathbb E\big[|Y_n|\mathbf 1_{\{|Y_n|>K\}}\big]=0.
$$
A standard theorem says:

**Vitali theorem.** If $X_n\xrightarrow{\mathbb P}X$ and $\{X_n\}$ is uniformly integrable, then $X\in L^1$ and
$$
X_n\xrightarrow{L^1}X.
$$

Useful sufficient conditions for uniform integrability are:

- domination: $|X_n|\le Y$ for all $n$ with $Y\in L^1$,
- an $L^{1+\delta}$ bound: $\sup_n \mathbb E|X_n|^{1+\delta}<\infty$ for some $\delta>0$.

#### Domination upgrades probability to \(L^p\)

If $X_n\xrightarrow{\mathbb P}X$ and $|X_n|\le Y$ almost surely for all $n$, where $Y\in L^p$, then in fact
$$
X_n\xrightarrow{L^p}X.
$$
A convenient proof uses the subsequence principle plus dominated convergence along almost surely convergent subsequences.

#### Summable error probabilities force almost sure convergence

If for each $\varepsilon>0$,
$$
\sum_n \mathbb P(|X_n-X|>\varepsilon)<\infty,
$$
then $X_n\to X$ almost surely.

A common sufficient condition is
$$
\sum_n \mathbb E|X_n-X|^p<\infty
\quad\Rightarrow\quad
X_n\to X \text{ almost surely},
$$
because Markov’s inequality gives summability of the tail probabilities.

---

## 8. Uniqueness and stability of limits

### Uniqueness up to almost sure equality

If $X_n\to X$ and $X_n\to Y$ in any of the modes

- almost surely,
- in probability,
- in $L^p$,

then $X=Y$ almost surely.

For convergence in probability, the proof is standard: for any $\varepsilon>0$,
$$
\mathbb P(|X-Y|>\varepsilon)
\le
\mathbb P(|X_n-X|>\varepsilon/2)+\mathbb P(|X_n-Y|>\varepsilon/2)\to 0.
$$
Hence $\mathbb P(|X-Y|>\varepsilon)=0$ for all $\varepsilon>0$, so $X=Y$ almost surely.

For almost sure convergence, uniqueness follows pointwise on the full-probability set where both limits exist. For $L^p$, uniqueness follows because $L^p$ convergence implies convergence in probability.

### Uniqueness of weak limits

If
$$
X_n\xrightarrow{d}X
\quad\text{and}\quad
X_n\xrightarrow{d}Y,
$$
then $X$ and $Y$ have the same distribution:
$$
\mathcal L(X)=\mathcal L(Y).
$$
Weak limits are therefore unique **in law**, not necessarily as random variables on a common space.

### Passing to subsequences

If $X_n\to X$ in any of the modes above, then every subsequence $X_{n_k}$ converges to the same limit in the same mode.

This is often used together with the subsequence principle: to prove convergence in probability, it is enough to show that every subsequence has an almost surely convergent subsubsequence with the right limit.

---

## 9. Standard theorems used all the time

### Continuous Mapping Theorem

Let $g:\mathbb R^m\to\mathbb R^k$ be measurable, and let $D_g$ be its discontinuity set.

**Weak version.** If $X_n\xrightarrow{d}X$ in $\mathbb R^m$ and
$$
\mathbb P(X\in D_g)=0,
$$
then
$$
g(X_n)\xrightarrow{d}g(X).
$$

There are analogous versions for almost sure convergence and convergence in probability:

- if $X_n\xrightarrow{\text{a.s.}}X$ and $g$ is continuous at $X(\omega)$ for almost every $\omega$, then
  $$
  g(X_n)\xrightarrow{\text{a.s.}}g(X);
  $$
- if $X_n\xrightarrow{\mathbb P}X$ and $g$ is continuous at $X$ almost surely, then
  $$
  g(X_n)\xrightarrow{\mathbb P}g(X).
  $$

### Slutsky’s theorem

A very common form is:

**Theorem (Slutsky).** If
$$
X_n\xrightarrow{d}X,
\qquad
Y_n\xrightarrow{\mathbb P}c,
$$
where $c$ is constant, then
$$
(X_n,Y_n)\xrightarrow{d}(X,c).
$$
Consequently,
$$
X_n+Y_n\xrightarrow{d}X+c,
\qquad
X_nY_n\xrightarrow{d}cX,
$$
and if $c\neq 0$,
$$
\frac{X_n}{Y_n}\xrightarrow{d}\frac{X}{c}.
$$

More generally, if $g$ is continuous at every point of the form $(x,c)$, then
$$
g(X_n,Y_n)\xrightarrow{d}g(X,c).
$$

This is really the combination of joint convergence with the Continuous Mapping Theorem.

### Delta method

Suppose $T_n$ is an estimator of $\theta$, and
$$
a_n(T_n-\theta)\xrightarrow{d}Z
$$
for some deterministic sequence $a_n\to\infty$.

If $g:\mathbb R\to\mathbb R$ is differentiable at $\theta$, then
$$
a_n\big(g(T_n)-g(\theta)\big)\xrightarrow{d}g'(\theta)\,Z.
$$

In several dimensions, if $T_n\in\mathbb R^d$, $g:\mathbb R^d\to\mathbb R^k$ is differentiable at $\theta$, and $Dg(\theta)$ is its Jacobian, then
$$
a_n\big(g(T_n)-g(\theta)\big)\xrightarrow{d}Dg(\theta)\,Z.
$$

The proof is a first-order Taylor expansion plus Slutsky.

### Subsequence principle

Already stated above:
$$
X_n\xrightarrow{\mathbb P}X
\quad\Longleftrightarrow\quad
\text{every subsequence has an a.s.-convergent subsubsequence}.
$$
This is one of the most useful structural facts about convergence in probability.

### Cramér–Wold device

For random vectors $X_n,X$ in $\mathbb R^d$,
$$
X_n\xrightarrow{d}X
\quad\Longleftrightarrow\quad
t^\top X_n\xrightarrow{d}t^\top X
\quad\text{for every }t\in\mathbb R^d.
$$
This is the standard way to prove weak convergence of random vectors.

### Skorokhod representation theorem

If $X_n\xrightarrow{d}X$ in a Polish space (for example, $\mathbb R^d$), then there exist random variables $Y_n,Y$ on a common probability space such that

- $\mathcal L(Y_n)=\mathcal L(X_n)$,
- $\mathcal L(Y)=\mathcal L(X)$,
- $Y_n\xrightarrow{\text{a.s.}}Y$.

This is extremely useful conceptually: weak convergence can be represented as almost sure convergence on another probability space.

But one must remember the caveat: **this new coupling is not the original one**. Skorokhod representation does **not** say that the original $X_n$ converge almost surely.

---

## 10. Tools for proving convergence

### Markov and Chebyshev inequalities

For $p>0$ and any random variable $Y$,
$$
\mathbb P(|Y|>\varepsilon)\le \frac{\mathbb E|Y|^p}{\varepsilon^p}.
$$
This is the quickest proof that $L^p$ convergence implies convergence in probability.

For $p=2$, one gets Chebyshev’s inequality:
$$
\mathbb P(|Y-\mathbb EY|>\varepsilon)\le \frac{\operatorname{Var}(Y)}{\varepsilon^2}.
$$
This is often how weak laws of large numbers are proved.

### Borel–Cantelli lemmas

If $\sum_n \mathbb P(A_n)<\infty$, then
$$
\mathbb P(A_n\ \text{i.o.})=0.
$$
If the $A_n$ are independent and $\sum_n \mathbb P(A_n)=\infty$, then
$$
\mathbb P(A_n\ \text{i.o.})=1.
$$

These lemmas are the basic bridge from probability bounds to almost sure statements.

### Dominated convergence and bounded convergence

If $Y_n\to Y$ almost surely and $|Y_n|\le Z$ with $Z\in L^1$, then
$$
\mathbb E[Y_n]\to \mathbb E[Y].
$$

A particularly useful corollary: if $X_n\to X$ almost surely and $|X_n|\le Y$ with $Y\in L^p$, then
$$
X_n\xrightarrow{L^p}X,
$$
because $|X_n-X|^p\le (2Y)^p$ and $|X_n-X|^p\to 0$ almost surely.

### Uniform integrability

Uniform integrability is the standard condition that allows expectation-type convergence to pass through convergence in probability.

As above, if $X_n\xrightarrow{\mathbb P}X$ and $\{X_n\}$ is uniformly integrable, then
$$
X_n\xrightarrow{L^1}X.
$$

This theorem explains why convergence of distributions alone is not enough to pass to expectations of unbounded functions.

### Tightness on \(\mathbb R\)

A sequence $\{X_n\}$ is **tight** if for every $\varepsilon>0$ there exists $M<\infty$ such that
$$
\sup_n \mathbb P(|X_n|>M)<\varepsilon.
$$
If $X_n\xrightarrow{d}X$, then $\{X_n\}$ is tight.

On $\mathbb R$, tightness is the main compactness condition for weak convergence. In practice, a common strategy is:

1. prove tightness,
2. identify all possible subsequential weak limits,
3. conclude that the whole sequence converges in distribution.

---

## 11. Algebra of convergent random variables

### Almost sure and probability convergence

If
$$
X_n\to X,\qquad Y_n\to Y
$$
almost surely, then
$$
(X_n,Y_n)\to (X,Y)
$$
almost surely, and therefore for every continuous $g:\mathbb R^2\to\mathbb R$,
$$
g(X_n,Y_n)\to g(X,Y)
$$
almost surely.

The same is true for convergence in probability: if
$$
X_n\xrightarrow{\mathbb P}X,\qquad Y_n\xrightarrow{\mathbb P}Y,
$$
then
$$
(X_n,Y_n)\xrightarrow{\mathbb P}(X,Y),
$$
because for the sup norm,
$$
\mathbb P\!\left(\max\{|X_n-X|,|Y_n-Y|\}>\varepsilon\right)
\le
\mathbb P(|X_n-X|>\varepsilon)+\mathbb P(|Y_n-Y|>\varepsilon)\to 0.
$$
Hence, by the Continuous Mapping Theorem for probability convergence,
$$
g(X_n,Y_n)\xrightarrow{\mathbb P}g(X,Y)
$$
for continuous $g$.

This immediately yields the usual algebra:
$$
X_n+Y_n\to X+Y,\qquad
X_nY_n\to XY,\qquad
|X_n|\to |X|,
$$
and
$$
\max(X_n,Y_n)\to \max(X,Y),\qquad
\min(X_n,Y_n)\to \min(X,Y)
$$
in the same mode (a.s. or in probability).

For quotients, one needs the denominator not to hit $0$ in the limit. For example, if
$$
X_n\xrightarrow{\mathbb P}X,\qquad Y_n\xrightarrow{\mathbb P}Y,\qquad \mathbb P(Y=0)=0,
$$
then
$$
\frac{X_n}{Y_n}\xrightarrow{\mathbb P}\frac{X}{Y}.
$$

### Convergence in distribution

For weak convergence, functions of several variables require **joint** convergence. If
$$
(X_n,Y_n)\xrightarrow{d}(X,Y),
$$
then for every continuous $g$,
$$
g(X_n,Y_n)\xrightarrow{d}g(X,Y).
$$

A common mistake is to think that separate convergence
$$
X_n\xrightarrow{d}X,\qquad Y_n\xrightarrow{d}Y
$$
automatically implies
$$
(X_n,Y_n)\xrightarrow{d}(X,Y).
$$
This is false in general.

Slutsky is the important special case where one component converges to a constant.

### \(L^p\) algebra

For $p\ge 1$,

- if $X_n\xrightarrow{L^p}X$ and $Y_n\xrightarrow{L^p}Y$, then
  $$
  X_n+Y_n\xrightarrow{L^p}X+Y
  $$
  by Minkowski’s inequality;

- if $g$ is Lipschitz, then
  $$
  \|g(X_n)-g(X)\|_p\le L\|X_n-X\|_p,
  $$
  so
  $$
  X_n\xrightarrow{L^p}X
  \quad\Longrightarrow\quad
  g(X_n)\xrightarrow{L^p}g(X).
  $$

Products require extra integrability. A standard useful fact is:

If
$$
X_n\xrightarrow{L^2}X,\qquad Y_n\xrightarrow{L^2}Y,
$$
then
$$
X_nY_n\xrightarrow{L^1}XY.
$$
Indeed,
$$
\|X_nY_n-XY\|_1
\le
\|X_n(Y_n-Y)\|_1+\|(X_n-X)Y\|_1
\le
\|X_n\|_2\|Y_n-Y\|_2+\|Y\|_2\|X_n-X\|_2,
$$
and the right-hand side tends to $0$ because $L^2$ convergence implies boundedness in $L^2$.

---

## 12. Asymptotic notation used in statistics

### Definitions

We write
$$
X_n=o_p(1)
$$
if
$$
X_n\xrightarrow{\mathbb P}0.
$$

We write
$$
X_n=O_p(1)
$$
if the sequence is bounded in probability: for every $\varepsilon>0$ there exist $M<\infty$ and $N$ such that
$$
\mathbb P(|X_n|>M)<\varepsilon
\qquad \text{for all }n\ge N.
$$
Equivalently, $\{X_n\}$ is tight.

More generally, for deterministic $a_n>0$,
$$
X_n=o_p(a_n)
\quad\Longleftrightarrow\quad
\frac{X_n}{a_n}\xrightarrow{\mathbb P}0,
$$
and
$$
X_n=O_p(a_n)
\quad\Longleftrightarrow\quad
\frac{X_n}{a_n}=O_p(1).
$$

### Basic algebra

The following rules are used constantly:

- $o_p(1)+o_p(1)=o_p(1)$,
- $O_p(1)+O_p(1)=O_p(1)$,
- $O_p(1)\,o_p(1)=o_p(1)$,
- if $X_n\xrightarrow{d}X$, then $X_n=O_p(1)$,
- if $X_n\xrightarrow{\mathbb P}c\neq 0$, then
  $$
  \frac1{X_n}\xrightarrow{\mathbb P}\frac1c.
  $$

### How Slutsky is used in practice

A typical argument is:
$$
T_n = Z_n + o_p(1),
\qquad
Z_n\xrightarrow{d}Z.
$$
Then
$$
T_n\xrightarrow{d}Z.
$$

Another common pattern is studentization:
if
$$
\sqrt n(\hat\theta_n-\theta)\xrightarrow{d}Z,
\qquad
\hat s_n\xrightarrow{\mathbb P}s>0,
$$
then
$$
\frac{\sqrt n(\hat\theta_n-\theta)}{\hat s_n}
\xrightarrow{d}
\frac{Z}{s}.
$$

---

## 13. Worked examples

### Example 1: independent rare events

Let $X_n\sim \mathrm{Bernoulli}(1/n)$ independently.

- For any $\varepsilon\in(0,1)$,
  $$
  \mathbb P(|X_n|>\varepsilon)=\frac1n\to 0,
  $$
  so $X_n\to 0$ in probability.

- For any $p>0$,
  $$
  \mathbb E|X_n|^p=\frac1n\to 0,
  $$
  so $X_n\to 0$ in $L^p$ for every $p>0$.

- But
  $$
  \sum_n \mathbb P(X_n=1)=\infty,
  $$
  and independence gives
  $$
  \mathbb P(X_n=1\ \text{i.o.})=1,
  $$
  so $X_n\not\to 0$ almost surely.

This single example shows that even very strong $L^p$ convergence need not imply almost sure convergence.

### Example 2: deterministic spikes

Let
$$
X_n=n\,\mathbf 1_{(0,1/n)} \quad \text{on }([0,1],\mathcal B,\lambda).
$$

- For each $\omega>0$, $X_n(\omega)=0$ eventually, so $X_n\to 0$ almost surely.
- Therefore $X_n\to 0$ in probability and in distribution.
- But
  $$
  \mathbb E|X_n|=1,
  $$
  so there is no $L^1$ convergence.
- Also
  $$
  \mathbb E[X_n]=1 \not\to 0=\mathbb E[0].
  $$

This is the standard warning that weak/probability convergence does not control expectations of unbounded functions.

### Example 3: same law, no probability convergence

Let $X$ be Rademacher and define $X_n=(-1)^nX$.

- Every $X_n$ has the same law as $X$, so $X_n\xrightarrow{d}X$.
- But $X_n$ has no probability limit.

This shows that convergence in distribution only concerns marginals, not the coupling across $n$.

### Example 4: \(L^1\) but not \(L^2\)

Let
$$
X_n=\sqrt n\,\mathbf 1_{(0,1/n)}.
$$
Then
$$
\mathbb E|X_n|=\frac1{\sqrt n}\to 0
$$
but
$$
\mathbb E(X_n^2)=1.
$$
Hence $X_n\to 0$ in $L^1$ but not in $L^2$.

### Example 5: law of large numbers and central limit theorem

Let $X_1,X_2,\dots$ be i.i.d. with mean $\mu$ and variance $\sigma^2<\infty$, and define
$$
\bar X_n = \frac1n\sum_{i=1}^n X_i.
$$

- **Weak law of large numbers:**
  $$
  \bar X_n\xrightarrow{\mathbb P}\mu.
  $$
- **Strong law of large numbers:** if $\mathbb E|X_1|<\infty$, then
  $$
  \bar X_n\xrightarrow{\text{a.s.}}\mu.
  $$
- **Central limit theorem:**
  $$
  \sqrt n(\bar X_n-\mu)\xrightarrow{d}N(0,\sigma^2).
  $$

If $S_n^2$ is the sample variance and $S_n^2\xrightarrow{\mathbb P}\sigma^2$, then by Slutsky,
$$
\frac{\sqrt n(\bar X_n-\mu)}{S_n}\xrightarrow{d}N(0,1).
$$

This is the standard asymptotic-statistics pipeline:
strong/weak law for consistency, then CLT plus Slutsky for asymptotic normality of standardized statistics.

---

## 14. Summary sheets

### 14.1 One table of all definitions

| Mode            | Notation                        | Definition                                                   | Same probability space needed? |
| --------------- | ------------------------------- | ------------------------------------------------------------ | ------------------------------ |
| Almost sure     | $X_n\xrightarrow{\text{a.s.}}X$ | $\mathbb P(\{\omega:X_n(\omega)\to X(\omega)\})=1$           | Yes                            |
| In probability  | $X_n\xrightarrow{\mathbb P}X$   | $\mathbb P(|X_n-X|>\varepsilon)\to 0$ for every $\varepsilon>0$ | Yes                            |
| In $L^p$        | $X_n\xrightarrow{L^p}X$         | $\mathbb E|X_n-X|^p\to 0$                                    | Yes                            |
| In distribution | $X_n\xrightarrow{d}X$           | $F_{X_n}(x)\to F_X(x)$ at every continuity point of $F_X$    | No                             |

### 14.2 One table of equivalent definitions

| Mode            | Equivalent formulations most worth remembering               |
| --------------- | ------------------------------------------------------------ |
| Almost sure     | $\mathbb P(|X_n-X|>\varepsilon\ \text{i.o.})=0$ for every $\varepsilon>0$; equivalently $\{X_n\not\to X\}=\bigcup_{m\ge1}\limsup_n\{|X_n-X|>1/m\}$ |
| In probability  | convergence in measure; metric convergence under $d(U,V)=\mathbb E[\min\{|U-V|,1\}]$; every subsequence has an a.s.-convergent subsubsequence |
| In $L^p$        | $\mathbb E|X_n-X|^p\to 0$; for $p\ge1$, $\|X_n-X\|_p\to 0$   |
| In distribution | bounded continuous test functions; Portmanteau theorem; characteristic functions (Lévy continuity theorem) |

### 14.3 One implication chart

| Implication                            | True? | Why / counterexample                                         |
| -------------------------------------- | ----- | ------------------------------------------------------------ |
| $L^q\Rightarrow L^p$ for $q\ge p$      | Yes   | Norm monotonicity on probability spaces                      |
| $L^p\Rightarrow \mathbb P$             | Yes   | Markov inequality                                            |
| a.s. $\Rightarrow \mathbb P$           | Yes   | Dominated convergence on indicators                          |
| $\mathbb P\Rightarrow d$               | Yes   | Bounded continuous test functions / subsequence principle    |
| Complete $\Rightarrow$ a.s.            | Yes   | Borel–Cantelli                                               |
| $d\Rightarrow \mathbb P$               | No    | $X_n=(-1)^nX$ with $X$ Rademacher                            |
| $\mathbb P\Rightarrow$ a.s.            | No    | Independent Bernoulli$(1/n)$                                 |
| a.s. $\Rightarrow L^1$                 | No    | $n\,\mathbf 1_{(0,1/n)}$                                     |
| $L^1\Rightarrow$ a.s.                  | No    | Independent Bernoulli$(1/n)$                                 |
| $L^1\Rightarrow L^2$                   | No    | $\sqrt n\,\mathbf 1_{(0,1/n)}$                               |
| $d\Rightarrow \mathbb P$ to a constant | Yes   | Weak convergence to a constant equals probability convergence |

