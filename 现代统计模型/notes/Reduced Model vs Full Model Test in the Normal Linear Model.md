# Reduced Model vs Full Model Test in the Normal Linear Model
## A Detailed Derivation via Orthogonal Decomposition and Residual Regression

## 1. Problem setup

Consider the normal linear model
$$
y = X\beta + \varepsilon, \qquad \varepsilon \sim N(0,\sigma^2 I_n),
$$
where:
- $y \in \mathbb R^n$,
- $X \in \mathbb R^{n\times p}$,
- $\beta \in \mathbb R^p$,
- $\sigma^2 > 0$.

We assume that the column space of the reduced model is contained in that of the full model:
$$
\mathcal C(X_0) \subset \mathcal C(X),
$$
where $X_0 \in \mathbb R^{n\times p_0}$ is the design matrix of the reduced model.

We want to test
$$
H_0: E(y)\in \mathcal C(X_0)
\qquad \text{vs} \qquad
H_A: E(y)\in \mathcal C(X)\setminus \mathcal C(X_0).
$$

A standard way to represent this is to write the full model as
$$
X = [X_0,\; X_1],
$$
so that
$$
y = X_0\beta_0 + X_1\beta_1 + \varepsilon.
$$
Then the hypothesis is informally “the additional part $X_1$ contributes nothing beyond $X_0$.” 
However, because the columns of $X_1$ may overlap with those of $X_0$, the correct object is not $X_1$ itself, but the part of $X_1$ orthogonal to $\mathcal C(X_0)$.

---

## 2. Why we do not literally use $\mathcal C(X)\setminus \mathcal C(X_0)$

The set
$$
\mathcal C(X)\setminus \mathcal C(X_0)
$$
is not a linear subspace, so it is not suitable for linear-model decomposition. 
What we actually want is the subspace of directions in $\mathcal C(X)$ that are orthogonal to $\mathcal C(X_0)$:
$$
\mathcal E := \mathcal C(X)\cap \mathcal C(X_0)^\perp.
$$

This is the genuinely “new” part of the full model beyond the reduced model.
The goal is to show that
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal E,
$$
and then reduce the hypothesis test to testing whether the mean has any component in $\mathcal E$.

---

## 3. Projection onto the reduced-model space

Let
$$
P_0 := P_{\mathcal C(X_0)}
$$
denote the orthogonal projection onto $\mathcal C(X_0)$. In matrix form,
$$
P_0 = X_0(X_0'X_0)^{-}X_0',
$$
where $(X_0'X_0)^{-}$ denotes any generalized inverse if needed.

Recall the basic properties of an orthogonal projection:
$$
P_0' = P_0, \qquad P_0^2 = P_0.
$$
Hence
$$
I-P_0
$$
is also symmetric and idempotent, and it is the orthogonal projection onto $\mathcal C(X_0)^\perp$.

Now define
$$
R_1 := (I-P_0)X_1.
$$
This is the residual of $X_1$ after projecting out its component in $\mathcal C(X_0)$.

---

## 4. The key decomposition of $\mathcal C(X)$

We now prove that
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C\big((I-P_0)X_1\big).
$$

### 4.1 First inclusion: $\mathcal C(X)\subset \mathcal C(X_0)+\mathcal C((I-P_0)X_1)$

Take any vector $v\in \mathcal C(X)$. Then there exist vectors $a,b$ such that
$$
v = X_0 a + X_1 b.
$$
Now decompose $X_1b$ into its projection onto $\mathcal C(X_0)$ and its orthogonal residual:
$$
X_1b = P_0X_1b + (I-P_0)X_1b.
$$
Since $P_0X_1b \in \mathcal C(X_0)$, we get
$$
v = \underbrace{(X_0a + P_0X_1b)}_{\in \mathcal C(X_0)}
+ \underbrace{(I-P_0)X_1b}_{\in \mathcal C((I-P_0)X_1)}.
$$
Therefore
$$
v \in \mathcal C(X_0)+\mathcal C((I-P_0)X_1).
$$
Since $v$ was arbitrary,
$$
\mathcal C(X)\subset \mathcal C(X_0)+\mathcal C((I-P_0)X_1).
$$

### 4.2 Reverse inclusion: $\mathcal C(X_0)+\mathcal C((I-P_0)X_1)\subset \mathcal C(X)$

First, clearly
$$
\mathcal C(X_0)\subset \mathcal C(X).
$$
Also,
$$
(I-P_0)X_1 = X_1 - P_0X_1.
$$
Here $X_1 \in \mathcal C(X)$ columnwise, and $P_0X_1 \in \mathcal C(X_0)\subset \mathcal C(X)$. Hence
$$
(I-P_0)X_1 \in \mathcal C(X)
$$
columnwise, so
$$
\mathcal C((I-P_0)X_1)\subset \mathcal C(X).
$$
Thus
$$
\mathcal C(X_0)+\mathcal C((I-P_0)X_1)\subset \mathcal C(X).
$$

Combining both inclusions,
$$
\mathcal C(X)=\mathcal C(X_0)+\mathcal C((I-P_0)X_1).
$$

### 4.3 Orthogonality

We now show the sum is orthogonal.
Take any $u\in \mathcal C(X_0)$ and any $w\in \mathcal C((I-P_0)X_1)$. Then there exist vectors $c,d$ such that
$$
u = X_0 c, \qquad w = (I-P_0)X_1 d.
$$
Then
$$
u'w = c'X_0'(I-P_0)X_1 d.
$$
But since every column of $X_0$ lies in $\mathcal C(X_0)$,
$$
P_0X_0 = X_0.
$$
Therefore
$$
X_0'(I-P_0)=X_0' - X_0'P_0 = X_0' - X_0' = 0.
$$
Hence
$$
u'w = 0.
$$
So
$$
\mathcal C(X_0)\perp \mathcal C((I-P_0)X_1).
$$

Therefore we have the orthogonal direct-sum decomposition
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C((I-P_0)X_1).
$$

---

## 5. Identification of the “new” subspace

By construction,
$$
\mathcal C((I-P_0)X_1)\subset \mathcal C(X_0)^\perp,
$$
because $(I-P_0)$ maps into $\mathcal C(X_0)^\perp$.

Also, as shown above,
$$
\mathcal C((I-P_0)X_1)\subset \mathcal C(X).
$$
Thus
$$
\mathcal C((I-P_0)X_1)\subset \mathcal C(X)\cap \mathcal C(X_0)^\perp.
$$

Conversely, let $v\in \mathcal C(X)\cap \mathcal C(X_0)^\perp$. Since $v\in\mathcal C(X)$, and we already proved
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C((I-P_0)X_1),
$$
there exist unique $u\in \mathcal C(X_0)$ and $w\in \mathcal C((I-P_0)X_1)$ such that
$$
v=u+w.
$$
Because $v\perp \mathcal C(X_0)$ and $w\perp \mathcal C(X_0)$, we must also have
$$
u = v-w \perp \mathcal C(X_0).
$$
But $u\in \mathcal C(X_0)$ as well, so $u=0$. Therefore $v=w\in \mathcal C((I-P_0)X_1)$.

Hence
$$
\mathcal C((I-P_0)X_1)=\mathcal C(X)\cap \mathcal C(X_0)^\perp.
$$

So the “extra” subspace is exactly
$$
\mathcal E = \mathcal C((I-P_0)X_1).
$$

---

## 6. Rewriting the full model using the orthogonal decomposition

Since the full-model mean space decomposes as
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal E,
\qquad \mathcal E=\mathcal C((I-P_0)X_1),
$$
every mean vector in $\mathcal C(X)$ can be written uniquely as
$$
\mu = X_0\alpha + E\gamma,
$$
where $E$ is any matrix whose column space is $\mathcal E$. A natural choice is
$$
E=(I-P_0)X_1.
$$
Hence the full model may be rewritten as
$$
y = X_0\alpha + E\gamma + \varepsilon.
$$

Under this parameterization:
- the reduced model is exactly the set of mean vectors with $\gamma=0$,
- the full model allows arbitrary $\gamma$.

Therefore
$$
H_0: E(y)\in \mathcal C(X_0)
\qquad \Longleftrightarrow \qquad
H_0:\gamma=0.
$$

This is the precise sense in which the original reduced-vs-full test becomes a test on the orthogonal complement.

---

## 7. Residualizing the response

Now apply $(I-P_0)$ to the full model:
$$
(I-P_0)y = (I-P_0)X_0\alpha + (I-P_0)E\gamma + (I-P_0)\varepsilon.
$$
Since $X_0$ lies in $\mathcal C(X_0)$,
$$
(I-P_0)X_0 = 0.
$$
Since $E\in \mathcal C(X_0)^\perp$, we have
$$
(I-P_0)E = E.
$$
Therefore
$$
(I-P_0)y = E\gamma + (I-P_0)\varepsilon.
$$

Define
$$
r_y := (I-P_0)y.
$$
Then
$$
r_y = E\gamma + (I-P_0)\varepsilon.
$$

This is the residualized model. It says:

> After removing the part explained by $X_0$, the remaining signal is exactly the component in the extra subspace $\mathcal E$.

Since $E=(I-P_0)X_1$, we may also write
$$
r_y = (I-P_0)X_1\gamma + (I-P_0)\varepsilon.
$$

Thus the original hypothesis test is equivalent to testing
$$
H_0:\gamma=0
$$
in the regression of the residualized response on the residualized extra regressors.

---

## 8. Distribution of the residualized error

Because
$$
\varepsilon \sim N(0,\sigma^2 I_n),
$$
and $(I-P_0)$ is a fixed symmetric idempotent matrix,
$$
(I-P_0)\varepsilon \sim N\big(0,\sigma^2(I-P_0)\big).
$$

So the residualized model is still Gaussian:
$$
r_y \sim N\big(E\gamma,\sigma^2(I-P_0)\big).
$$
On the subspace $\mathcal C(X_0)^\perp$, the covariance acts like $\sigma^2 I$, so ordinary least squares on the residualized regression is valid.

---

## 9. Projection matrices for the reduced and full models

Let
$$
P_X := P_{\mathcal C(X)}
$$
denote the projection onto the full-model space.

Since
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(E)
$$
with orthogonal summands, the projection onto $\mathcal C(X)$ is the sum of the two orthogonal projections:
$$
P_X = P_0 + P_E,
$$
where
$$
P_E := P_{\mathcal C(E)} = P_{\mathcal C((I-P_0)X_1)}.
$$
Thus
$$
P_X - P_0 = P_E.
$$

This identity is fundamental. It says that the difference between the full-model projection and the reduced-model projection is exactly the projection onto the new subspace contributed by the full model.

---

## 10. Sum-of-squares decomposition

The residual sum of squares under the reduced model is
$$
SSE_R = y'(I-P_0)y.
$$
The residual sum of squares under the full model is
$$
SSE_F = y'(I-P_X)y.
$$

Using $P_X=P_0+P_E$, we get
$$
I-P_X = I-P_0-P_E.
$$
Hence
$$
SSE_F = y'(I-P_0-P_E)y.
$$
Subtracting,
$$
SSE_R - SSE_F
= y'(I-P_0)y - y'(I-P_0-P_E)y
= y'P_E y.
$$
Since $P_E(I-P_0)=P_E$ and $(I-P_0)P_E=P_E$,
$$
y'P_E y = y'(I-P_0)P_E(I-P_0)y.
$$
Therefore
$$
SSE_R - SSE_F = r_y'P_E r_y.
$$

This is exactly the regression sum of squares when regressing $r_y$ on $E$.

So the extra sum of squares due to adding $X_1$ beyond $X_0$ is the sum of squares explained by the residualized extra regressors in the residualized response.

---

## 11. Rank and degrees of freedom

Let
$$
r_0 := \operatorname{rank}(X_0), \qquad r := \operatorname{rank}(X).
$$
Then, because
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(E),
$$
we have
$$
\operatorname{rank}(E)=r-r_0.
$$
Set
$$
q := r-r_0.
$$
This is the number of additional linearly independent directions introduced by the full model beyond the reduced model.

The error degrees of freedom under the full model are
$$
n-r.
$$

---

## 12. Distribution of the quadratic forms under $H_0$

Under $H_0$, the mean vector satisfies
$$
E(y)\in \mathcal C(X_0).
$$
Therefore the component in $\mathcal C(E)$ is zero, so
$$
P_E E(y)=0.
$$
Hence under $H_0$,
$$
y'P_E y = \varepsilon'P_E\varepsilon.
$$
Since $P_E$ is symmetric idempotent of rank $q$,
$$
\frac{y'P_E y}{\sigma^2} \sim \chi^2_q.
$$

Likewise, under the full model residual,
$$
SSE_F = y'(I-P_X)y.
$$
Because $I-P_X$ is symmetric idempotent of rank $n-r$ and annihilates the mean under the full model,
$$
\frac{SSE_F}{\sigma^2} \sim \chi^2_{n-r}.
$$

Moreover, because
$$
P_E(I-P_X)=0,
$$
the two projection matrices project onto orthogonal subspaces, so the corresponding quadratic forms are independent:
$$
y'P_E y \;\perp\; y'(I-P_X)y.
$$

Therefore under $H_0$,
$$
\frac{(SSE_R-SSE_F)/q}{SSE_F/(n-r)}
=
\frac{(y'P_E y)/q}{y'(I-P_X)y/(n-r)}
\sim F_{q,n-r}.
$$

This is the partial $F$-test for comparing the reduced and full models.

---

## 13. The multi-step regression interpretation

The above derivation justifies the following three-step procedure.

### Step 1: Regress $X_1$ on $X_0$

Compute the residual matrix
$$
R_1=(I-P_0)X_1.
$$
This removes the part of $X_1$ already explained by $X_0$.

### Step 2: Regress $y$ on $X_0$

Compute the residual vector
$$
r_y=(I-P_0)y.
$$
This removes the part of $y$ already explained by $X_0$.

### Step 3: Regress $r_y$ on $R_1$

Fit
$$
r_y = R_1\gamma + \text{error}.
$$
Then test
$$
H_0:\gamma=0.
$$

This regression has exactly the same numerator sum of squares, denominator sum of squares, and therefore the same $F$-statistic as the original reduced-vs-full comparison.

So the “multi-step regression method” is not a heuristic; it is an exact reformulation of the reduced-vs-full test based on the orthogonal decomposition of the mean space.

---

## 14. Relation to the Frisch–Waugh–Lovell idea

This derivation is a geometric version of the Frisch–Waugh–Lovell theorem.

The theorem says that the effect of $X_1$ after controlling for $X_0$ can be obtained by:
1. removing the effect of $X_0$ from $X_1$,
2. removing the effect of $X_0$ from $y$,
3. regressing the second residual on the first.

In our setting, the theorem is not only about coefficient estimation; it also gives the exact testing decomposition underlying the partial $F$-test.

---

## 15. A simple scalar example

Suppose the full model is
$$
y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \varepsilon,
$$
and the reduced model is
$$
y = \beta_0 + \beta_1 x_1 + \varepsilon.
$$

Then:
- $X_0=[\mathbf 1,\; x_1]$,
- $X_1=[x_2]$.

The new part of the full model is not $x_2$ itself, but
$$
r_{x_2} = (I-P_0)x_2,
$$
the residual from regressing $x_2$ on $(1,x_1)$.

Similarly,
$$
r_y=(I-P_0)y
$$
is the residual from regressing $y$ on $(1,x_1)$.

Then the original hypothesis
$$
H_0:\beta_2=0 \text{ (after controlling for intercept and }x_1)
$$
is equivalent to testing
$$
H_0:\gamma=0
$$
in
$$
r_y = \gamma r_{x_2} + \text{error}.
$$

This makes the interpretation especially clear:

> $x_2$ is significant only if the part of $x_2$ that is unrelated to $(1,x_1)$ explains the part of $y$ that is unrelated to $(1,x_1)$.

---

## 16. Final summary

The reduced-vs-full model test can be derived rigorously as follows:

1. Start from
$$
\mathcal C(X_0)\subset \mathcal C(X).
$$

2. Define the projection onto the reduced-model space:
$$
P_0 = P_{\mathcal C(X_0)}.
$$

3. Extract the orthogonal complement contributed by the full model:
$$
\mathcal E = \mathcal C((I-P_0)X_1).
$$

4. Prove the orthogonal decomposition:
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal E.
$$

5. Rewrite the model as
$$
y = X_0\alpha + E\gamma + \varepsilon,
\qquad E=(I-P_0)X_1.
$$

6. Reduce the hypothesis to
$$
H_0:\gamma=0.
$$

7. Residualize both response and added regressors:
$$
r_y=(I-P_0)y,\qquad R_1=(I-P_0)X_1.
$$

8. Regress
$$
r_y \sim R_1
$$
and test whether the coefficient is zero.

9. The corresponding extra sum of squares is
$$
SSE_R-SSE_F = y'(P_X-P_0)y = y'P_E y.
$$

10. Under $H_0$,
$$
F=\frac{(SSE_R-SSE_F)/q}{SSE_F/(n-r)}\sim F_{q,n-r}.
$$

So the multi-step regression method is an exact and geometrically transparent reformulation of the partial $F$-test: it isolates the component of the full model that is genuinely new relative to the reduced model and tests whether the mean has any component in that direction.



# Reduced vs Full Model Test via Pure Matrix Methods

## 1. Setup

Consider the normal linear model
$$
y = X\beta + \varepsilon, \qquad \varepsilon \sim N(0,\sigma^2 I_n).
$$

Assume
$$
\mathcal C(X_0) \subset \mathcal C(X),
$$
and we want to test
$$
H_0: E(y)\in \mathcal C(X_0)
\qquad \text{vs} \qquad
H_A: E(y)\in \mathcal C(X)\setminus \mathcal C(X_0).
$$

Write the full design matrix as
$$
X = [X_0,\; X_1].
$$

Let
$$
P_0 := P_{\mathcal C(X_0)}, \qquad P := P_{\mathcal C(X)}
$$
be the orthogonal projection matrices onto $\mathcal C(X_0)$ and $\mathcal C(X)$, respectively.

---

## 2. Goal of the matrix derivation

We want to derive, using only matrix arguments, that:

1. the "new part" of the full model beyond the reduced model is
$$
\mathcal C\big((I-P_0)X_1\big),
$$

2. the full-model space decomposes as
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C\big((I-P_0)X_1\big),
$$

3. the reduced-vs-full test is equivalent to a regression of
$$
(I-P_0)y
$$
on
$$
(I-P_0)X_1,
$$

4. the extra sum of squares is
$$
SSE_R-SSE_F = y'(P-P_0)y,
$$

5. under $H_0$, the test statistic
$$
F = \frac{y'(P-P_0)y / q}{y'(I-P)y / (n-r)}
$$
has an $F_{q,n-r}$ distribution, where
$$
q = \operatorname{rank}(X)-\operatorname{rank}(X_0), \qquad r=\operatorname{rank}(X).
$$

---

## 3. Basic matrix facts

We repeatedly use the following facts.

### 3.1 Orthogonal projections

If $P_M$ is the orthogonal projector onto a subspace $M$, then
$$
P_M' = P_M, \qquad P_M^2 = P_M.
$$

Also,
- $\mathcal C(P_M)=M$,
- $\mathcal C(I-P_M)=M^\perp$.

### 3.2 Characterization by action on subspaces

If $M_1 \subset M_2$, then
$$
P_{M_2}P_{M_1} = P_{M_1}, \qquad P_{M_1}P_{M_2}=P_{M_1}.
$$

Since $\mathcal C(X_0)\subset \mathcal C(X)$, we have
$$
PP_0 = P_0P = P_0.
$$

### 3.3 Orthogonal direct sum and projections

If
$$
M = M_1 \oplus M_2
$$
with $M_1 \perp M_2$, then
$$
P_M = P_{M_1}+P_{M_2}.
$$

---

## 4. Constructing the new part of the full model

Define
$$
R_1 := (I-P_0)X_1.
$$

We claim that
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(R_1).
$$

### 4.1 Show $\mathcal C(X)\subset \mathcal C(X_0)+\mathcal C(R_1)$

Take any vector $v\in\mathcal C(X)$. Since $X=[X_0,X_1]$, there exist vectors $a,b$ such that
$$
v = X_0 a + X_1 b.
$$

Now insert the decomposition
$$
X_1 b = P_0X_1 b + (I-P_0)X_1 b.
$$
Then
$$
v = \big(X_0 a + P_0X_1 b\big) + (I-P_0)X_1 b.
$$

Here,
$$
X_0 a + P_0X_1 b \in \mathcal C(X_0),
$$
and
$$
(I-P_0)X_1 b \in \mathcal C(R_1).
$$

Hence
$$
v\in \mathcal C(X_0)+\mathcal C(R_1).
$$

Therefore
$$
\mathcal C(X)\subset \mathcal C(X_0)+\mathcal C(R_1).
$$

### 4.2 Show the reverse inclusion

Clearly
$$
\mathcal C(X_0)\subset \mathcal C(X).
$$

Also,
$$
R_1 = (I-P_0)X_1 = X_1 - P_0X_1.
$$
Since both $X_1$ and $P_0X_1$ have columns in $\mathcal C(X)$, each column of $R_1$ lies in $\mathcal C(X)$. Hence
$$
\mathcal C(R_1)\subset \mathcal C(X).
$$

Therefore
$$
\mathcal C(X_0)+\mathcal C(R_1)\subset \mathcal C(X).
$$

Combining the two inclusions gives
$$
\mathcal C(X)=\mathcal C(X_0)+\mathcal C(R_1).
$$

### 4.3 Show orthogonality

Take $u\in\mathcal C(X_0)$ and $w\in\mathcal C(R_1)$. Then
$$
u=X_0c, \qquad w=(I-P_0)X_1 d
$$
for some vectors $c,d$. Thus
$$
u'w = c'X_0'(I-P_0)X_1 d.
$$

Since the columns of $X_0$ lie in $\mathcal C(X_0)$,
$$
P_0X_0 = X_0,
$$
so
$$
X_0'(I-P_0)=0.
$$
Hence
$$
u'w=0.
$$

Therefore
$$
\mathcal C(X_0)\perp \mathcal C(R_1),
$$
and so
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(R_1).
$$

This is the key geometric decomposition, obtained purely by matrix manipulations.

---

## 5. The projection identity $P-P_0 = P_{R_1}$

Let
$$
P_{R_1} := P_{\mathcal C(R_1)}.
$$

Since we just proved
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(R_1),
$$
and the sum is orthogonal, the projector onto $\mathcal C(X)$ is
$$
P = P_0 + P_{R_1}.
$$

Therefore
$$
P - P_0 = P_{R_1}.
$$

This identity is central: it says that the additional projection gained by moving from the reduced model to the full model is exactly the projection onto the residualized regressor space.

---

## 6. Residual regression emerges algebraically

The reduced-model residual vector is
$$
r_y := (I-P_0)y.
$$

Because $P_{R_1}$ projects onto a subspace of $\mathcal C(X_0)^\perp$, we have
$$
P_{R_1}(I-P_0)=P_{R_1}, \qquad (I-P_0)P_{R_1}=P_{R_1}.
$$

Thus
$$
P_{R_1}y = P_{R_1}(I-P_0)y = P_{R_1}r_y.
$$

So the part of $y$ additionally explained by the full model beyond the reduced model depends only on the reduced-model residual $r_y$, and it is obtained by projecting $r_y$ onto $\mathcal C(R_1)$.

Since
$$
R_1=(I-P_0)X_1,
$$
this is exactly the regression of $(I-P_0)y$ on $(I-P_0)X_1$.

Hence the multi-step regression procedure is a purely matrix consequence of the projection decomposition:
1. residualize $X_1$:
$$
R_1=(I-P_0)X_1,
$$
2. residualize $y$:
$$
r_y=(I-P_0)y,
$$
3. regress $r_y$ on $R_1$.

No heuristic argument is needed.

---

## 7. Sum-of-squares decomposition by matrices alone

The reduced-model residual sum of squares is
$$
SSE_R = y'(I-P_0)y.
$$

The full-model residual sum of squares is
$$
SSE_F = y'(I-P)y.
$$

Subtract:
$$
SSE_R - SSE_F
= y'(I-P_0)y - y'(I-P)y
= y'(P-P_0)y.
$$

Using $P-P_0=P_{R_1}$, we obtain
$$
SSE_R - SSE_F = y'P_{R_1}y.
$$

Now use the fact that $P_{R_1}(I-P_0)=P_{R_1}$:
$$
y'P_{R_1}y
= y'(I-P_0)P_{R_1}(I-P_0)y
= r_y'P_{R_1}r_y.
$$

Thus
$$
SSE_R-SSE_F = r_y'P_{R_1}r_y.
$$

This is exactly the regression sum of squares from regressing $r_y$ on $R_1$.

So the extra sum of squares is not merely analogous to the residual regression; it is literally the same quadratic form.

---

## 8. Algebraic form of the hypothesis

The null hypothesis is
$$
H_0: E(y)\in\mathcal C(X_0).
$$

Since $\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(R_1)$, every vector in $\mathcal C(X)$ can be uniquely written as
$$
\mu = \mu_0 + \mu_1,
\qquad
\mu_0\in\mathcal C(X_0), \quad \mu_1\in\mathcal C(R_1).
$$

Therefore
$$
H_0 \iff \mu_1 = 0.
$$

Equivalently, since $P_{R_1}$ extracts precisely the $\mathcal C(R_1)$-component,
$$
H_0 \iff P_{R_1}E(y)=0.
$$

This is the matrix formulation of “the additional directions contribute nothing.”

---

## 9. Ranks and degrees of freedom

Let
$$
r_0 := \operatorname{rank}(X_0), \qquad r:=\operatorname{rank}(X).
$$

Because
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(R_1),
$$
we have
$$
\operatorname{rank}(R_1)=r-r_0.
$$

Define
$$
q := r-r_0.
$$

Then
- $q$ is the number of extra linearly independent directions in the full model,
- $n-r$ is the residual degrees of freedom under the full model.

---

## 10. Distribution of the numerator quadratic form

Under $H_0$, we have $E(y)\in\mathcal C(X_0)$. Since $\mathcal C(R_1)\perp\mathcal C(X_0)$,
$$
P_{R_1}E(y)=0.
$$
Thus under $H_0$,
$$
y'P_{R_1}y = \varepsilon'P_{R_1}\varepsilon.
$$

Now $P_{R_1}$ is symmetric and idempotent with rank $q$. A standard quadratic-form result gives
$$
\frac{\varepsilon'P_{R_1}\varepsilon}{\sigma^2}\sim \chi_q^2.
$$

Therefore
$$
\frac{SSE_R-SSE_F}{\sigma^2}
=
\frac{y'P_{R_1}y}{\sigma^2}
\sim \chi_q^2.
$$

---

## 11. Distribution of the denominator quadratic form

Similarly,
$$
SSE_F = y'(I-P)y.
$$

Because $I-P$ is symmetric and idempotent with rank $n-r$, and because $(I-P)E(y)=0$ whenever $E(y)\in\mathcal C(X)$, we have
$$
\frac{SSE_F}{\sigma^2}\sim \chi_{n-r}^2.
$$

---

## 12. Independence of numerator and denominator

To obtain the $F$ distribution, we also need independence.

Observe that
$$
P_{R_1}(I-P)=0,
$$
because $\mathcal C(R_1)\subset \mathcal C(X)$, while $I-P$ projects onto $\mathcal C(X)^\perp$.

Since the two projection matrices are symmetric, idempotent, and project onto orthogonal subspaces, the corresponding Gaussian quadratic forms are independent. Hence
$$
y'P_{R_1}y
\quad \text{and} \quad
y'(I-P)y
$$
are independent under $H_0$.

---

## 13. The $F$ statistic

Therefore under $H_0$,
$$
F
=
\frac{(SSE_R-SSE_F)/q}{SSE_F/(n-r)}
=
\frac{y'(P-P_0)y/q}{y'(I-P)y/(n-r)}
\sim F_{q,n-r}.
$$

Using $P-P_0=P_{R_1}$, we may also write
$$
F
=
\frac{y'P_{R_1}y/q}{y'(I-P)y/(n-r)}.
$$

This is the standard reduced-vs-full partial $F$ test, derived entirely from matrix identities.

---

## 14. Why the multi-step regression is exactly equivalent

The equivalence is now immediate:

- the extra part of the full fit is the projection onto $\mathcal C(R_1)$,
- this projection acts only on the reduced residual $r_y=(I-P_0)y$,
- therefore the gain from the full model over the reduced model is exactly the fit obtained by regressing $r_y$ on $R_1=(I-P_0)X_1$.

So the following two procedures are mathematically identical:

### Direct reduced-vs-full comparison
Compare
$$
SSE_R = y'(I-P_0)y
\quad \text{and} \quad
SSE_F = y'(I-P)y.
$$

### Multi-step regression comparison
Regress
$$
r_y=(I-P_0)y
$$
on
$$
R_1=(I-P_0)X_1
$$
and use its regression sum of squares.

They produce the same numerator, the same denominator, and hence the same $F$ statistic.

---

## 15. Final conclusion

Yes—the entire method can be derived by pure matrix techniques.

The essential chain is:

1. Define the reduced-model projector $P_0$.
2. Residualize the added regressors:
$$
R_1=(I-P_0)X_1.
$$
3. Prove the orthogonal decomposition
$$
\mathcal C(X)=\mathcal C(X_0)\oplus \mathcal C(R_1).
$$
4. Deduce the projection identity
$$
P-P_0=P_{R_1}.
$$
5. Rewrite the extra sum of squares as
$$
SSE_R-SSE_F = y'(P-P_0)y = r_y'P_{R_1}r_y.
$$
6. Use quadratic-form theory to show
$$
\frac{SSE_R-SSE_F}{\sigma^2}\sim \chi_q^2,
\qquad
\frac{SSE_F}{\sigma^2}\sim \chi_{n-r}^2,
$$
and that these are independent under $H_0$.
7. Conclude
$$
F=\frac{(SSE_R-SSE_F)/q}{SSE_F/(n-r)}\sim F_{q,n-r}.
$$

So the multi-step regression method is not just an intuitive construction; it is the exact matrix-algebraic form of the reduced-vs-full test.